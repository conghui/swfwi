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
#include "essfwiframework.h"
#include "shotdata-reader.h"
#include "sfutil.h"
#include "sum.h"
#include "Matrix.h"
#include "calgainmatrix.h"
#include "preserved-alpha.h"

static const int N = 2;
static const int enkf_update_every_essfwi_iter = 1;

//static void initVelSet(const fwi_config_t &config, std::vector<float *> &velSet, int N, int modelSize) {
//  TRACE() << "add perturbation to initial velocity";
//  const std::string fn = "velPerturb_" + to_str(N) + ".bin";
//  std::ifstream ifs(fn.c_str());
//  if (!ifs) {
//    ERROR() << "cannot open file: " << fn;
//    exit(EXIT_FAILURE);
//  }
//
//  std::vector<float> tmp(modelSize);
//  std::vector<float> velOrig(config.init_vel, config.init_vel + modelSize);
//  std::transform(velOrig.begin(), velOrig.end(), velOrig.begin(),
//                 boost::bind(velRecover<float>, _1, config.fwi_dim.dx, config.dt));
//
//  for (std::vector<float *>::iterator it = velSet.begin(); it != velSet.end(); ++it) {
//    float *ret = *it;
//    ifs.read(reinterpret_cast<char *>(&tmp[0]), modelSize * sizeof(tmp[0]));
//    std::transform(tmp.begin(), tmp.end(), velOrig.begin(), ret, std::plus<float>());
//    std::transform(ret, ret + modelSize, ret, boost::bind(velTrans<float>, _1, config.fwi_dim.dx, config.dt));
//  }
//  ifs.close();
//}
template <typename T>
std::string to_str(T val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}

template <typename E>
float velRecover(float vel, float dx, float dt) {
  float t = dx * dx / (dt * dt * vel);
  return std::sqrt(t);
}

template <typename E>
float velTrans(float vel, float dx, float dt) {
  return dx * dx / (dt * dt * vel * vel);
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

void enkf_analyze(Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, std::vector<float *> &velSet, int modelSize, int iter) {
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

  TRACE() << "go through ENKF to update velocity set";
  enkf_analyze(fmMethod, wlt, dobs, velSet, modelSize, 1);

  float *vel = createAMean(velSet, modelSize);
  writeVelocity("meamvel0.rsf", vel, exvel.nx, exvel.nz, dx, dt);
  finalizeAMean(vel);

  std::vector<Damp4t10d *> fms(N);
  std::vector<EssFwiFrameworkOld *> essfwis(N);
  for (size_t i = 0; i < essfwis.size(); i++) {
    fms[i] = new Damp4t10d(fmMethod);
    fms[i]->bindVelocity(*veldb[i]);
    essfwis[i] = new EssFwiFrameworkOld(*fms[i], wlt, dobs);
  }

  TRACE() << "iterate the remaining iteration";
  for (int iter = 1; iter <= params.niter; iter++) {
    TRACE() << "FWI for each velocity";
    DEBUG() << "\n\n\n\n\n\n\n";

    for (int ivel = 0; ivel < N; ivel++) {
      essfwis[ivel]->epoch(iter, ivel);
//      epoch(config, velSet[i], NULL, curr_gradient, update_direction, iter, i + 1, N, gradUpdator);
    }
    TRACE() << "enkf analyze and update velocity";
    if (iter % enkf_update_every_essfwi_iter == 0) {
      enkf_analyze(fmMethod, wlt, dobs, velSet, modelSize, iter + 1);

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
