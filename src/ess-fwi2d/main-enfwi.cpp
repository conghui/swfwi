extern "C" {
#include <rsf.h>
}

#include <cstdlib>

#include "logger.h"
#include "essfwi-params.h"
#include "shot-position.h"
#include "damp4t10d.h"
#include "sf-velocity-reader.h"
#include "ricker-wavelet.h"
#include "essfwiframework.h"
#include "shotdata-reader.h"
#include "sum.h"

static const int N = 20;

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

static std::vector<float *> createVelocitySet(int modelSize, int num) {
  std::vector<float *> ret(num);

  for (size_t i = 0; i < ret.size(); i++) {
    ret[i] = (float *)malloc(modelSize * sizeof * ret[i]);
    std::fill(ret[i], ret[i] + modelSize, 0);
  }

  return ret;
}

static void initVelSetRandom(const Velocity &vel, std::vector<float *> &velSet) {
  int modelSize = vel.nx * vel.nz;
  std::vector<float> tmp(modelSize);
  const std::vector<float> &velOrig = vel.dat;
//  std::transform(velOrig.begin(), velOrig.end(), velOrig.begin(), boost::bind(velRecover<float>, _1, config.fwi_dim.dx, config.dt));

  for (std::vector<float *>::iterator it = velSet.begin(); it != velSet.end(); ++it) {
    float *ret = *it;
    float noise = (rand() % 10 + 1) / 10f;
    std::fill(tmp.begin(), tmp.end(), noise);
    std::transform(tmp.begin(), tmp.end(), velOrig.begin(), ret, std::plus<float>());
//    std::transform(ret, ret + modelSize, ret, boost::bind(velTrans<float>, _1, config.fwi_dim.dx, config.dt));
  }
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

void enkf_analyze(Damp4t10d &fm, std::vector<float *> &velSet, int modelSize, int iter) {
  std::vector<float *> &A = velSet; /// velSet <==> matA
  int N = velSet.size();
  float *AMean = createAMean(A, modelSize);
  DEBUG() << "sum of AMean: " << std::accumulate(AMean, AMean + modelSize, 0.0f);

//  Matrix A_Perturb(N, modelSize);
//  initAPerturb(A_Perturb, A, AMean, modelSize);
//  DEBUG() << "sum of A_Perturb: " << getSum(A_Perturb);
//
//  std::vector<float> resdSet(N);
//
//  Matrix gainMatrix = calGainMatrix(config, velSet, resdSet, AMean, modelSize, iter);
//
//  Matrix t5(N, modelSize);
//  alpha_A_B_plus_beta_C(1, A_Perturb, gainMatrix, 0, t5);
//  DEBUG() << "sum of t5: " << getSum(t5);
//
//  TRACE() << "add the update back to velocity model";
//  for (size_t i = 0; i < velSet.size(); i++) {
//    float *vel = velSet[i];
//
//    DEBUG() << format("before vel recovery, velset[%2d/%d], min: %f, max: %f") % i % velSet.size() %
//            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));
//
//    TRACE() << "transform velocity to original";
//    std::transform(vel, vel + modelSize, vel, boost::bind(velRecover<float>, _1, config.fwi_dim.dx, config.dt));
//
//    DEBUG() << format("vel recoverty,       velset[%2d/%d], min: %f, max: %f") % i % velSet.size() %
//            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));
//
//    TRACE() << "add value calculated from ENKF to velocity";
//    Matrix::value_type *pu = t5.getData() + i * t5.getNumRow();
//    std::transform(vel, vel + modelSize, pu, vel, std::plus<Matrix::value_type>());
//
//    TRACE() << "sum velset " << i << ": " << std::accumulate(vel, vel + modelSize, 0.0f);
//    DEBUG() << format("after plus ENKF,     velset[%2d/%d], min: %f, max: %f\n") % i % velSet.size() %
//            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));
//
//    std::transform(vel, vel + modelSize, vel, boost::bind(velTrans<float>, _1, config.fwi_dim.dx, config.dt));
//  }

  //    TRACE() << "updating lamdaset";
  //    Matrix Lambda_Perturb(N, 2);
  //    initLambdaPerturb(lambdaSet, Lambda_Perturb);
  //
  //    Matrix t6(N, 2);
  //    alpha_A_B_plus_beta_C(1, Lambda_Perturb, gainMatrix, 0, t6);
  //
  //    A_plus_B(lambdaSet, t6, lambdaSet);

//  TRACE() << "updating ratioset";
//  Matrix ratio_Perturb(N, 2);
//  initRatioPerturb(ratioSet, ratio_Perturb);
//
//  Matrix t6(N, 2);
//  alpha_A_B_plus_beta_C(1, ratio_Perturb, gainMatrix, 0, t6);
//
//  A_plus_B(ratioSet, t6, ratioSet);
//
//  TRACE() << "updating lamdaset";
//  const dim3d_t &dim = config.fwi_dim;
//  for (int i = 0; i < N; i++) {
//    ReguFactor fac(velSet[i], dim.nx, dim.nz, config.lambdaX, config.lambdaZ);
//    Matrix::value_type *pLambda = lambdaSet.getData() + i * lambdaSet.getNumRow();
//    Matrix::value_type *pRatio  = ratioSet.getData()  + i * ratioSet.getNumRow();
//    pLambda[0] = pRatio[0] * resdSet[i] / fac.getWx2();
//    pLambda[1] = pRatio[1] * resdSet[i] / fac.getWz2();
//    DEBUG() << "PRATIO[0] " << pRatio[0];
//    DEBUG() << "PRESD[i] " << resdSet[i];
//    DEBUG() << "PLAMBDA[0] " << pLambda[0];
//    DEBUG() << "";
//  }

  //    TRACE() << "Another forwarding to test D-HA and D-HAa";
  //    if (true) {
  //      for (size_t i = 0; i < velSet.size(); i++) {
  //        DEBUG() << format("calculate HA on velocity %2d/%d") % (i + 1) % velSet.size();
  //        TRACE() << "making encoded shot, both sources and receivers";
  //
  //        shot_t shot;
  //        TRACE() << "TESTING: every column use same code, codes[0]";
  //        if (diffColDiffCodes) {
  //          TRACE() << "different codes for different column";
  //          make_encoded_shot(config, config.shots_config.n_shots, &shot, codes[i]);
  //        } else {
  //          TRACE() << "same codes for different column";
  //          make_encoded_shot(config, config.shots_config.n_shots, &shot, codes[0]);
  //        }
  //
  //        TRACE() << "expand velocity and shot";
  //        dim3d_t dim = config.fwi_dim;
  //        expand_velocity_model(dim, &velSet[i], config.freesurface, 1);
  //        expand_shot(&shot, dim, 1);
  //        expand_dim(&dim, config.freesurface, 1);
  //
  //        TRACE() << "H operate on A, and store data in HOnA";
  //        forward_modeling_dispatch(dim, shot, velSet[i], config.forward_modeling_boundary, config.check_step, config.use_gpu);
  //
  //        const Matrix::value_type *pd = D.getData() + i * numDataSamples;
  //        const Matrix::value_type *pha = HOnA.getData() + i * numDataSamples;
  //        const float *phaa = shot.shot_output.traces[0].data;
  //        float varDHA = 0;
  //        float varDHAa = 0;
  //        for (int j = 0; j < numDataSamples; j++) {
  //          varDHA += (pd[j] - pha[j]) * (pd[j] - pha[j]);
  //          varDHAa += (pd[j] - phaa[j]) * (pd[j] - phaa[j]);
  //        }
  //        DEBUG() << format("norm2(D-HA) = %f, norm2(D-HAa): %f") % varDHA % varDHAa;
  //
  //        TRACE() << "shrink the dimensions back to original";
  //        shrink_dim(&dim, config.freesurface, 1);
  //        shrink_velocity_model(dim, &velSet[i], config.freesurface, 1);
  //
  //        free_shot(&shot);
  //      } /// end of for
  //    } /// end of if
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
  std::vector<float *> velSet = createVelocitySet(modelSize, N);
  initVelSetRandom(exvel, velSet); /// init the perturbation within the program later

  TRACE() << "go through ENKF to update velocity set";
//  enkf_analyze(config, velSet, lambdaSet, ratioSet, modelSize, 1);


  EssFwiFramework essfwi(fmMethod, wlt, dobs);

  for (int iter = 0; iter < params.niter; iter++) {
    essfwi.epoch(iter);
    essfwi.writeVel(params.vupdates);
  } /// end of iteration

  sf_close();

  return 0;
}
