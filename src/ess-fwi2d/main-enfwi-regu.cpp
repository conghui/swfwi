extern "C" {
#include <rsf.h>
}

#include <cstdlib>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <string>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "logger.h"
#include "essfwi-params.h"
#include "shot-position.h"
#include "damp4t10d.h"
#include "sf-velocity-reader.h"
#include "ricker-wavelet.h"
#include "essfwiregu.h"
#include "shotdata-reader.h"
#include "sfutil.h"
#include "sum.h"
#include "Matrix.h"
#include "calgainmatrix.h"
#include "preserved-alpha.h"
#include "ReguFactor.h"
#include "random-code.h"
#include "encoder.h"
#include "common.h"

static const int N = 2;
static const int enkf_update_every_essfwi_iter = 1;
static const float initLambdaRatio = 0.5;

template <typename T>
std::string to_str(T val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}


template <typename T>
float variance(const T *A_Begin, const T *A_End, const T *B_Begin) {
  float v = 0;

  for (; A_Begin != A_End; ++A_Begin, ++B_Begin) {
    v += (*A_Begin - *B_Begin) * (*A_Begin - *B_Begin);
  }

  return v;
}


template <typename T>
T variance(const std::vector<T> &A, const std::vector<T> &B) {
  const T *A_begin = &A[0];
  const T *A_end   = &A[A.size() - 1] + 1;
  const T *B_begin = &B[0];
  return variance(A_begin, A_end, B_begin);
}

static float *createAMean(const std::vector<float *> &velSet, int modelSize) {
  float *ret = (float *)malloc(modelSize * sizeof * ret);

  for (int i = 0; i < modelSize; i++) {
    float sum = 0.0;
    for (size_t j = 0; j < velSet.size(); j++) {
      sum += velSet[j][i];
    }

    ret[i] = sum / velSet.size();
  }

  return ret;
}

static void initAPerturb(Matrix &matAPerturb, const std::vector<float *> &velSet,
    const float *AMean, int modelSize) {
  assert((size_t)matAPerturb.getNumCol() == velSet.size()); // the matrix is row-majored
  assert(matAPerturb.getNumRow() == modelSize);

  for (int i = 0; i < matAPerturb.getNumCol(); i++) {
    Matrix::value_type *p = matAPerturb.getData() + (i * matAPerturb.getNumRow());
    std::transform(velSet[i], velSet[i] + modelSize, AMean, p, std::minus<Matrix::value_type>());
  }
}
static void writeVelocity(const std::string &fn, const float *cur_vel, int nx, int nz, float dx, float dt) {
  int tmp_size = nx * nz;
  std::vector<float> tmp_vel(tmp_size);
  std::transform(cur_vel, cur_vel + tmp_size, tmp_vel.begin(), boost::bind(velRecover<float>, _1, dx, dt));

  sfFloatWrite2d(fn.c_str(), &tmp_vel[0], nz, nx);
}

void initRatioPerturb(const Matrix &ratioSet, Matrix &ratioPerturb) {
  initGamma(ratioSet, ratioPerturb);
}

void enkf_analyze(Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, std::vector<float *> &velSet,
    Matrix &lambdaSet, Matrix &ratioSet, int modelSize, int iter) {
  float dt = fm.getdt();
  float dx = fm.getdx();

  std::vector<float *> &A = velSet; /// velSet <==> matA
  int N = velSet.size();
  float *AMean = createAMean(A, modelSize);
  DEBUG() << "sum of AMean: " << std::accumulate(AMean, AMean + modelSize, 0.0f);

  Matrix A_Perturb(N, modelSize);
  initAPerturb(A_Perturb, A, AMean, modelSize);
  DEBUG() << "sum of A_Perturb: " << getSum(A_Perturb);

  std::vector<float> resdSet(N);

  Matrix gainMatrix = calGainMatrix(fm, wlt, dobs, velSet, resdSet, AMean, modelSize, iter);

  Matrix t5(N, modelSize);
  alpha_A_B_plus_beta_C(1, A_Perturb, gainMatrix, 0, t5);
  DEBUG() << "sum of t5: " << getSum(t5);

  TRACE() << "add the update back to velocity model";
  for (size_t i = 0; i < velSet.size(); i++) {
    float *vel = velSet[i];

    DEBUG() << format("before vel recovery, velset[%2d/%d], min: %f, max: %f") % i % velSet.size() %
            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

    TRACE() << "transform velocity to original";
    std::transform(vel, vel + modelSize, vel, boost::bind(velRecover<float>, _1, dx, dt));

    DEBUG() << format("vel recoverty,       velset[%2d/%d], min: %f, max: %f") % i % velSet.size() %
            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

    TRACE() << "add value calculated from ENKF to velocity";
    Matrix::value_type *pu = t5.getData() + i * t5.getNumRow();
    std::transform(vel, vel + modelSize, pu, vel, std::plus<Matrix::value_type>());

    TRACE() << "sum velset " << i << ": " << std::accumulate(vel, vel + modelSize, 0.0f);
    DEBUG() << format("after plus ENKF,     velset[%2d/%d], min: %f, max: %f\n") % i % velSet.size() %
            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

    std::transform(vel, vel + modelSize, vel, boost::bind(velTrans<float>, _1, dx, dt));
  }

  TRACE() << "updating ratioset";
  Matrix ratio_Perturb(N, 2);
  initRatioPerturb(ratioSet, ratio_Perturb);

  Matrix t6(N, 2);
  alpha_A_B_plus_beta_C(1, ratio_Perturb, gainMatrix, 0, t6);

  A_plus_B(ratioSet, t6, ratioSet);

  int nx = fm.getnx();
  int nz = fm.getnz();
  float lambdaX = 0;
  float lambdaZ = 0;
  TRACE() << "updating lamdaset";
  for (int i = 0; i < N; i++) {
    ReguFactor fac(velSet[i], nx, nz, lambdaX, lambdaZ);
    Matrix::value_type *pLambda = lambdaSet.getData() + i * lambdaSet.getNumRow();
    Matrix::value_type *pRatio  = ratioSet.getData()  + i * ratioSet.getNumRow();
    pLambda[0] = pRatio[0] * resdSet[i] / fac.getWx2();
    pLambda[1] = pRatio[1] * resdSet[i] / fac.getWz2();
    DEBUG() << "PRATIO[0] " << pRatio[0];
    DEBUG() << "PRESD[i] " << resdSet[i];
    DEBUG() << "PLAMBDA[0] " << pLambda[0];
    DEBUG() << "";
  }

} /// end of function

std::vector<Velocity *> createVelDB(const Velocity &vel, int N, float dx, float dt) {
  std::vector<Velocity *> veldb(N);   /// here is all the velocity samples resides, others are pointer to this

  int modelSize = vel.nx * vel.nz;

  TRACE() << "add perturbation to initial velocity";
//  const std::string fn = "velPerturb_" + to_str(N) + ".bin";
  const std::string fn = "velPerturb_20.bin";
  std::ifstream ifs(fn.c_str());
  if (!ifs) {
    ERROR() << "cannot open file: " << fn;
    exit(EXIT_FAILURE);
  }

  std::vector<float> tmp(modelSize);
  std::vector<float> velOrig = vel.dat;
  std::transform(velOrig.begin(), velOrig.end(), velOrig.begin(), boost::bind(velRecover<float>, _1, dx, dt));
//  sfFloatWrite2d("initvel.rsf", &velOrig[0], vel.nz, vel.nx);

  for (int iv = 0; iv < N; iv++) {
    std::vector<float> ret(modelSize);
    ifs.read(reinterpret_cast<char *>(&tmp[0]), modelSize * sizeof(tmp[0]));
    std::transform(tmp.begin(), tmp.end(), velOrig.begin(), ret.begin(), std::plus<float>());
    std::transform(ret.begin(), ret.end(), ret.begin(), boost::bind(velTrans<float>, _1, dx, dt));

    veldb[iv] = new Velocity(ret, vel.nx, vel.nz);
  }

  ifs.close();

  return veldb;
}

std::vector<float *> generateVelSet(std::vector<Velocity *> &veldb) {
  std::vector<float *> velSet(veldb.size());
  for (size_t iv = 0; iv < veldb.size(); iv++) {
    velSet[iv] = &veldb[iv]->dat[0];
  }

  return velSet;
}

static void forwardModeling(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    std::vector<float> &dobs /* output (fast: ng, slow: nt) */)
{
    int nx = fmMethod.getnx();
    int nz = fmMethod.getnz();
    int ns = fmMethod.getns();
    int ng = fmMethod.getng();
    int nt = fmMethod.getnt();

    std::vector<float> p0(nz * nx, 0);
    std::vector<float> p1(nz * nx, 0);

    for(int it=0; it<nt; it++) {

      fmMethod.addEncodedSource(&p1[0], &encSrc[it * ns]);

      fmMethod.stepForward(&p0[0], &p1[0]);

      fmMethod.recordSeis(&dobs[it*ng], &p0[0]);

      std::swap(p1, p0);

    }

}

void initLambdaSet(const Damp4t10d &fm, const std::vector<float> &wlt,
    const std::vector<float> &dobs, std::vector<float *> &velSet, Matrix &lambdaSet,
                   const Matrix &ratioSet, int enkf_samples) {
  DEBUG() << "initial lambda ratio is: " << ratioSet.getData()[0];

  int ng = fm.getng();
  int nt = fm.getnt();
  int ns = fm.getns();
  int nx = fm.getnx();
  int nz = fm.getnz();
  int modelSize = nx * nz;
  int diffColDiffCodes = 0;

  int N = velSet.size();
  std::vector<std::vector<int> > codes(N);
  assert(N == enkf_samples);

  int numDataSamples = ng * nt;
  std::vector<float> obsData(numDataSamples);
  std::vector<float> synData(numDataSamples);

  for (int i = 0; i < enkf_samples; i++) {
    DEBUG() << format("init Lambda on velocity %2d/%d") % (i + 1) % velSet.size();

    TRACE() << "making encoded shot, both sources and receivers";
    codes[i] = RandomCode::genPlus1Minus1(ns);
    std::copy(codes[i].begin(), codes[i].end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

    Encoder encoder(codes[0]);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);


    TRACE() << "save encoded data";
    std::vector<float> trans(encobs.size());
    matrix_transpose(&encobs[0], &trans[0], ng, nt);
    const float *p= &trans[0];
    std::copy(p, p + numDataSamples, obsData.begin());

    TRACE() << "expand velocity and shot";

    TRACE() << "forward modeling";
    std::vector<float> dcal(encobs.size(), 0);

    ///// this decrease the performance ///
    Damp4t10d newfm = fm;
    Velocity curvel(std::vector<float>(velSet[i], velSet[i] + modelSize), fm.getnx(), fm.getnz());
    newfm.bindVelocity(curvel);
    forwardModeling(newfm, encsrc, dcal);
    matrix_transpose(&dcal[0], &trans[0], ng, nt);

    std::copy(p, p + numDataSamples, synData.begin());

    TRACE() << "calculate the data residule";
    float resd = variance(obsData, synData);

    TRACE() << "create regularization factor";
    ReguFactor reguFac(velSet[i], nx, nz);
    float Wx2 = reguFac.getWx2();
    float Wz2 = reguFac.getWz2();

    TRACE() << "update LambdaSet";
    assert(lambdaSet.getNumCol() == enkf_samples);
    Matrix::value_type *plamda = lambdaSet.getData() + i * lambdaSet.getNumRow();
    const Matrix::value_type *pratio = ratioSet.getData() + i * ratioSet.getNumRow();
    plamda[0] = pratio[0] * resd / Wx2; /// lambdaX
    plamda[1] = pratio[1] * resd / Wz2; /// lamdaZ
    DEBUG() << "pratio[0] " << pratio[0];
    DEBUG() << "plambda[0] " << plamda[0];

  }
} /// end of function


int main(int argc, char *argv[]) {

  /* initialize Madagascar */
  sf_init(argc, argv);

  Logger::instance().init("enfwi");

  EssFwiParams &params = EssFwiParams::instance();

  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float fm = params.fm;
  float dx = params.dx;

  // set random seed
  const int seed = 10;
  srand(seed);

  WARNING() << "For testing, we only set N to " << N;

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);
  Damp4t10d fmMethod(allSrcPos, allGeoPos, dt, dx, fm, nb, nt);

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);
  Velocity exvel = fmMethod.expandDomain(v0);
  fmMethod.bindVelocity(exvel);


  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  std::vector<float> dobs(ns * nt * ng);     /* all observed data */
  ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);

  int modelSize = fmMethod.getnx() * fmMethod.getnz();
  TRACE() << "init velocity set";
  DEBUG() << "model size " << modelSize;

  //////// use my functions
  std::vector<Velocity *> veldb = createVelDB(exvel, N, dx, dt);
  std::vector<float *> velSet = generateVelSet(veldb);
  ////////

  TRACE() << "Init ratio (mu) Set";
  Matrix ratioSet(N, 2);  /// 0 for muX, 1 for muZ
  std::fill(ratioSet.getData(), ratioSet.getData() + ratioSet.size(), initLambdaRatio);

  TRACE() << "init lambda set";
  Matrix lambdaSet(N, 2); /// 0 for lambdaX, 1 for lambdaZ
  initLambdaSet(fmMethod, wlt, dobs, velSet, lambdaSet, ratioSet, N);

  TRACE() << "go through ENKF to update velocity set";
  enkf_analyze(fmMethod, wlt, dobs, velSet, lambdaSet, ratioSet, modelSize, 1);

  float *vel = createAMean(velSet, modelSize);
  writeVelocity("meamvel0.rsf", vel, exvel.nx, exvel.nz, dx, dt);
  finalizeAMean(vel);

  std::vector<Damp4t10d *> fms(N);
  std::vector<EssFwiRegu *> essfwis(N);
  for (size_t i = 0; i < essfwis.size(); i++) {
    fms[i] = new Damp4t10d(fmMethod);
    fms[i]->bindVelocity(*veldb[i]);
    essfwis[i] = new EssFwiRegu(*fms[i], wlt, dobs);
  }

  TRACE() << "iterate the remaining iteration";
  for (int iter = 0; iter <= params.niter; iter++) {
    TRACE() << "FWI for each velocity";
    DEBUG() << "\n\n\n\n\n\n\n";

    for (int ivel = 0; ivel < N; ivel++) {
      const Matrix::value_type *lambda = lambdaSet.getData() + ivel * lambdaSet.getNumRow();
      float lambdaX = lambda[0];
      float lambdaZ = lambda[1];
      essfwis[ivel]->epoch(iter, ivel, lambdaX, lambdaZ);
//      epoch(config, velSet[i], NULL, curr_gradient, update_direction, iter, i + 1, N, gradUpdator);
    }
    TRACE() << "enkf analyze and update velocity";
    if (iter % enkf_update_every_essfwi_iter == 0) {
      enkf_analyze(fmMethod, wlt, dobs, velSet, lambdaSet, ratioSet, modelSize, iter + 1);

//      TRACE() << "assign the average of all stored alpha to each sample";
//      float *p = &PreservedAlpha::instance().getAlpha()[0];
//      float alphaAvg = std::accumulate(p, p + N, 0.0f) / N;
//      std::fill(p, p + N, alphaAvg);
    }

    float *vel = createAMean(velSet, modelSize);
    float l1norm, l2norm;
//    slownessL1L2Norm(config.accurate_vel, vel, config, l1norm, l2norm);
//    INFO() << format("%4d/%d iter, slowness l1norm: %g, slowness l2norm: %g") % iter % params.niter % l1norm % l2norm;

//    if (iter % 10 == 0) {
      TRACE() << "write the mean of velocity set";
      char buf[256];
      sprintf(buf, "vel%d.rsf", iter);
      writeVelocity(buf, vel, exvel.nx, exvel.nz, dx, dt);
//    }

    finalizeAMean(vel);
  }

  sf_close();

  return 0;
}
