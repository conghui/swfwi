/*
 * essfwiframework.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */


extern "C"
{
#include <rsf.h>
}

#include <time.h>
#include <cmath>

#include <omp.h>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cstdlib>
#include <functional>
#include <vector>
#include <set>

#include "logger.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "sum.h"
#include "sf-velocity-reader.h"
#include "shotdata-reader.h"
#include "random-code.h"
#include "encoder.h"
#include "velocity.h"
#include "sfutil.h"
#include "parabola-vertex.h"
#include "essfwiframework.h"

#include "aux.h"

namespace {

void updateGrad(float *pre_gradient, const float *cur_gradient, float *update_direction,
                           int model_size, int iter) {
  if (iter == 0) {
    std::copy(cur_gradient, cur_gradient + model_size, update_direction);
    std::copy(cur_gradient, cur_gradient + model_size, pre_gradient);
  } else {
    float beta = 0.0f;
    float a = 0.0f;
    float b = 0.0f;
    float c = 0.0f;
    int   i = 0;
    for (i = 0; i < model_size; i ++) {
      a += (cur_gradient[i] * cur_gradient[i]);
      b += (cur_gradient[i] * pre_gradient[i]);
      c += (pre_gradient[i] * pre_gradient[i]);
    }

    beta = (a - b) / c;

    if (beta < 0.0f) {
      beta = 0.0f;
    }

    for (i = 0; i < model_size; i ++) {
      update_direction[i] = cur_gradient[i] + beta * update_direction[i];
    }

    TRACE() << "Save current gradient to pre_gradient for the next iteration's computation";
    std::copy(cur_gradient, cur_gradient + model_size, pre_gradient);
  }
}


void second_order_virtual_source_forth_accuracy(float *vsrc, int num) {
  float *tmp_vsrc = (float *)malloc(num * sizeof(float));
  memcpy(tmp_vsrc, vsrc, num * sizeof(float));
  int i = 0;
  for (i = 0; i < num; i ++) {
    if ( i <= 1) {
      vsrc[i] = 0.0f;
      continue;
    }

    if ( (num - 1) == i || (num - 2) == i) {
      vsrc[i] = 0.0f;
      continue;
    }

    vsrc[i] = -1. / 12 * tmp_vsrc[i - 2] + 4. / 3 * tmp_vsrc[i - 1] -
              2.5 * tmp_vsrc[i] + 4. / 3 * tmp_vsrc[i + 1] - 1. / 12 * tmp_vsrc[i + 2];
  }

  free(tmp_vsrc);
}

void transVsrc(std::vector<float> &vsrc, int nt, int ng) {
  std::vector<float> trans(nt * ng);
  matrix_transpose(&vsrc[0], &trans[0], ng, nt);
  for (int ig = 0; ig < ng; ig++) {
    second_order_virtual_source_forth_accuracy(&trans[ig * nt], nt);
  }

  matrix_transpose(&trans[0], &vsrc[0], nt, ng);
}

void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale) {
  for (int i = 0; i < model_size; i ++) {
    image[i] -= src_wave[i] * vsrc_wave[i] * scale;
  }

}

void calgradient(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt)
{
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nz * nx, 0);
  std::vector<float> sp1(nz * nx, 0);
  std::vector<float> gp0(nz * nx, 0);
  std::vector<float> gp1(nz * nx, 0);


  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&sp1[0], &encSrc[it * ns], allSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);
    fmMethod.writeBndry(&bndr[0], &sp0[0], it);
  }

  for(int it = nt - 1; it >= 0 ; it--) {
    fmMethod.readBndry(&bndr[0], &sp0[0], it);
    std::swap(sp0, sp1);
    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    fmMethod.subEncodedSource(&sp0[0], &encSrc[it * ns]);

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0]);
    std::swap(gp1, gp0);

    if (dt * it > 0.4) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), 1.0);
    } else if (dt * it > 0.3) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), (dt * it - 0.3) / 0.1);
    } else {
      break;
    }
 }
}

void calgradient2_store(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt)
{
  const int check_step = 50;

  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nz * nx, 0);
  std::vector<float> sp1(nz * nx, 0);
  std::vector<float> gp0(nz * nx, 0);
  std::vector<float> gp1(nz * nx, 0);


  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&sp1[0], &encSrc[it * ns], allSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);
//    fmMethod.writeBndry(&bndr[0], &sp0[0], it);

    if ((it > 0) && (it != (nt - 1)) && !(it % check_step)) {
      char check_file_name1[64];
      char check_file_name2[64];
      const char *checkPointDir = std::getenv("CHECKPOINTDIR");
      sprintf(check_file_name1, "%s/check_time_%d_1.su", checkPointDir, it);
      sprintf(check_file_name2, "%s/check_time_%d_2.su", checkPointDir, it);
      writeBin(std::string(check_file_name1), &sp0[0], sp0.size() * sizeof(float));
      writeBin(std::string(check_file_name2), &sp1[0], sp1.size() * sizeof(float));
    }

  }

  char check_file_name1[64];
  char check_file_name2[64];
  const char *checkPointDir = std::getenv("CHECKPOINTDIR");
  sprintf(check_file_name1, "%s/check_time_last_1.su", checkPointDir);
  sprintf(check_file_name2, "%s/check_time_last_2.su", checkPointDir);
  writeBin(std::string(check_file_name1), &sp0[0], sp0.size() * sizeof(float));
  writeBin(std::string(check_file_name2), &sp1[0], sp1.size() * sizeof(float));

  for(int it = nt - 1; it >= 0 ; it--) {
//    fmMethod.addSource(&p1[0], &encSrc[it * ns], allSrcPos);
//    fmMethod.stepForward(&p0[0], &p1[0]);
//
    if (it  ==  nt - 1) {
      //Load last two time_step wave field
      char check_file_name1[64];
      char check_file_name2[64];
      const char *checkPointDir = std::getenv("CHECKPOINTDIR");
      sprintf(check_file_name1, "%s/check_time_last_1.su", checkPointDir);
      sprintf(check_file_name2, "%s/check_time_last_2.su", checkPointDir);
      readBin(std::string(check_file_name1), &sp1[0], sp1.size() * sizeof(float));
      readBin(std::string(check_file_name2), &sp0[0], sp0.size() * sizeof(float));
    }  else if ((check_step > 0) && !(it % check_step) && (it != 0)) {
      char check_file_name1[64];
      char check_file_name2[64];
      const char *checkPointDir = std::getenv("CHECKPOINTDIR");
      sprintf(check_file_name1, "%s/check_time_%d_1.su", checkPointDir, it);
      sprintf(check_file_name2, "%s/check_time_%d_2.su", checkPointDir, it);
      readBin(std::string(check_file_name1), &sp1[0], sp1.size() * sizeof(float));
      readBin(std::string(check_file_name2), &sp0[0], sp0.size() * sizeof(float));
//      printf("reading %s and %s\n", check_file_name1, check_file_name2);
    }

//    printf("it %d, check_step: %d\n", it, check_step);

//    {
//      char buf[256];
//      sprintf(buf, "sp1aaa%d.rsf", it);
//      sfFloatWrite2d(buf, &sp1[0], nzpad, nxpad);
//
//      sprintf(buf, "sp0aaa%d.rsf", it);
//      sfFloatWrite2d(buf, &sp0[0], nzpad, nxpad);
//    }

    fmMethod.stepBackward(&sp0[0], &sp1[0]);
    std::swap(sp1, sp0);

//    {
//      char buf[256];
//      sprintf(buf, "back%d.rsf", it);
//      sfFloatWrite2d(buf, &sp1[0], nzpad, nxpad);
//
//      sprintf(buf, "sp0back%d.rsf", it);
//      sfFloatWrite2d(buf, &sp0[0], nzpad, nxpad);
//    }

    fmMethod.subEncodedSource(&sp0[0], &encSrc[it * ns]);
//    fmMethod.subSource(&sp0[0], &encSrc[it * ns], all);
//    {
//      char buf[256];
//      sprintf(buf, "subw%d.rsf", it);
//      sfFloatWrite2d(buf, &sp0[0], nzpad, nxpad);
//    }

    /**
     * forward propagate receviers
     */
    fmMethod.addSource(&gp1[0], &vsrc[it * ng], allGeoPos);
//    {
//      char buf[256];
//      sprintf(buf, "vsrc%d.rsf", it);
//      sfFloatWrite1d(buf, &vsrc[it * ng], ng);
//    }
//    {
//      char buf[256];
//      sprintf(buf, "gp1aftadd%d.rsf", it);
//      sfFloatWrite2d(buf, &gp1[0], nzpad, nxpad);
//
//      sprintf(buf, "gp0aftadd%d.rsf", it);
//      sfFloatWrite2d(buf, &gp0[0], nzpad, nxpad);
//
//      sprintf(buf, "velaftadd%d.rsf", it);
//      sfFloatWrite2d(buf, &fmMethod.getVelocity().dat[0], nzpad, nxpad);
//    }

    fmMethod.stepForward(&gp0[0], &gp1[0]);
//    {
//      char buf[256];
//      sprintf(buf, "gp0aftfm%d.rsf", it);
//      sfFloatWrite2d(buf, &gp0[0], nzpad, nxpad);
//    }

    std::swap(gp1, gp0);

//    char buf[256];
//    sprintf(buf, "sfield%d.rsf", it);
//    sfFloatWrite2d(buf, &sp1[0], nzpad, nxpad);
//
//    sprintf(buf, "vfield%d.rsf", it);
//    sfFloatWrite2d(buf, &gp1[0], nzpad, nxpad);


    if (dt * it > 0.4) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), 1.0);
    } else if (dt * it > 0.3) {
      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), (dt * it - 0.3) / 0.1);
    } else {
      break;
    }

//    sprintf(buf, "img%d.rsf", it);
//    sfFloatWrite2d(buf, &g0[0], nzpad, nxpad);
//    if (it == 1999) exit(0);
 }

//  for(int it = nt - 1; it >= 0 ; it--) {
//    fmMethod.readBndry(&bndr[0], &sp0[0], it);
//    std::swap(sp0, sp1);
//    fmMethod.stepBackward(&sp0[0], &sp1[0]);
//    fmMethod.subEncodedSource(&sp0[0], &encSrc[it * ns]);
//
//    /**
//     * forward propagate receviers
//     */
//    fmMethod.addSource(&gp1[0], &vsrc[it * ng], allGeoPos);
//    fmMethod.stepForward(&gp0[0], &gp1[0]);
//    std::swap(gp1, gp0);
//
//    if (dt * it > 0.4) {
//      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), 1.0);
//    } else if (dt * it > 0.3) {
//      cross_correlation(&sp0[0], &gp0[0], &g0[0], g0.size(), (dt * it - 0.3) / 0.1);
//    } else {
//      break;
//    }
// }
}

} /// end of namespace


EssFwiFramework::EssFwiFramework(Damp4t10d &method, const UpdateSteplenOp &updateSteplenOp,
    const UpdateVelOp &_updateVelOp,
    const std::vector<float> &_wlt, const std::vector<float> &_dobs) :
    fmMethod(method), updateStenlelOp(updateSteplenOp), updateVelOp(_updateVelOp), wlt(_wlt), dobs(_dobs),
    ns(method.getns()), ng(method.getng()), nt(method.getnt()),
    nx(method.getnx()), nz(method.getnz()), dx(method.getdx()), dt(method.getdt()), objval(0)
{
  g0.resize(nx*nz, 0);
  updateDirection.resize(nx*nz, 0);
}

void EssFwiFramework::epoch(int iter) {
  // create random codes
  const std::vector<int> encodes = RandomCode::genPlus1Minus1(ns);

  Encoder encoder(encodes);
  std::vector<float> encsrc  = encoder.encodeSource(wlt);
  std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);

  std::vector<float> dcal(nt * ng, 0);
  fmMethod.EssForwardModeling(encsrc, dcal);
  fmMethod.removeDirectArrival(&encobs[0]);
  fmMethod.removeDirectArrival(&dcal[0]);

  std::vector<float> vsrc(nt * ng, 0);
  vectorMinus(encobs, dcal, vsrc);
  float obj1 = cal_objective(&vsrc[0], vsrc.size());
  DEBUG() << format("obj: %e") % obj1;

  transVsrc(vsrc, nt, ng);

  std::vector<float> g1(nx * nz, 0);
  calgradient(fmMethod, encsrc, vsrc, g1, nt, dt);

  DEBUG() << format("grad %.20f") % sum(g1);

  fmMethod.scaleGradient(&g1[0]);
  fmMethod.maskGradient(&g1[0]);

  updateGrad(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);

  updateStenlelOp.bindEncSrcObs(encsrc, encobs);
  float steplen;
  updateStenlelOp.calsteplen(updateDirection, obj1, iter, steplen, objval);

  Velocity &exvel = fmMethod.getVelocity();
  updateVelOp.update(exvel, exvel, updateDirection, steplen);

  fmMethod.refillBoundary(&exvel.dat[0]);
}

void EssFwiFramework::writeVel(sf_file file) const {
  fmMethod.sfWriteVel(fmMethod.getVelocity().dat, file);
}

float EssFwiFramework::getobjval() const {
  return objval;
}
