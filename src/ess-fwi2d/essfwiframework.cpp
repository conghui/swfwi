/*
 * essfwiframework.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#include "essfwiframework.h"



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

#include <boost/timer/timer.hpp>
#include "logger.h"
#include "essfwi-params.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "cycle-swap.h"
#include "sum.h"
#include "sf-velocity-reader.h"
#include "shotdata-reader.h"
#include "random-code.h"
#include "encoder.h"
#include "velocity.h"
#include "sfutil.h"
#include "aux.h"
#include "preserved-alpha.h"
#include "parabola-vertex.h"

std::vector<float> EssFwiFramework::g0;          /// gradient in previous step
std::vector<float> EssFwiFramework::updateDirection;

namespace {

static const int max_iter_select_alpha3 = 5;
static const float vmax = 5500;
static const float vmin = 1500;
static const float maxdv = 200;
typedef std::pair<float, float> ParaPoint;

bool parabolicLessComp(const ParaPoint &a, const ParaPoint &b) {
  return a.second - b.second < 1e-10;
}



void gradUpdator(float *, float *curr_grad, float *update_direction, int size) {
  std::copy(curr_grad, curr_grad + size, update_direction);
}

void prevCurrCorrDirection(float *pre_gradient, const float *cur_gradient, float *update_direction,
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


void forwardModeling(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    std::vector<float> &dobs, /* output (fast: ng, slow: nt) */
    int nt)
{
    int nxpad = fmMethod.getVelocity().nx;
    int nzpad = fmMethod.getVelocity().nz;
    int ns = fmMethod.getns();
    int ng = fmMethod.getng();

    boost::timer::cpu_timer timer;

    std::vector<float> p0(nzpad * nxpad, 0);
    std::vector<float> p1(nzpad * nxpad, 0);

    for(int it=0; it<nt; it++) {

      fmMethod.addEncodedSource(&p1[0], &encSrc[it * ns]);

      fmMethod.stepForward(&p0[0], &p1[0]);

      fmMethod.recordSeis(&dobs[it*ng], &p0[0]);

      std::swap(p1, p0);

    }

}

void vectorMinus(const std::vector<float> &dobs, const std::vector<float> &dcal, std::vector<float> &vsrc) {
  std::transform(dobs.begin(), dobs.end(), dcal.begin(), vsrc.begin(), std::minus<float>());
}

void second_order_virtual_source_forth_accuracy(float *vsrc, int num, float dt) {
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

void transVsrc(std::vector<float> &vsrc, int nt, int ng, float dt) {
  std::vector<float> trans(nt * ng);
  matrix_transpose(&vsrc[0], &trans[0], ng, nt);
  for (int ig = 0; ig < ng; ig++) {
    second_order_virtual_source_forth_accuracy(&trans[ig * nt], nt, dt);
  }

  //sfFloatWrite2d("vsrc.rsf", &trans[0], nt, ng);

  matrix_transpose(&trans[0], &vsrc[0], nt, ng);
}

static void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale) {
  int i = 0;
  for (i = 0; i < model_size; i ++) {
    image[i] -= src_wave[i] * vsrc_wave[i] * scale;
  }

}

void hello2(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt)
{
  const int check_step = 50;

  int nxpad = fmMethod.getnx();
  int nzpad = fmMethod.getnz();
  int ns = fmMethod.getns();
  int ng = fmMethod.getng();
  const ShotPosition &allGeoPos = fmMethod.getAllGeoPos();
  const ShotPosition &allSrcPos = fmMethod.getAllSrcPos();

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> sp0(nzpad * nxpad, 0);
  std::vector<float> sp1(nzpad * nxpad, 0);
  std::vector<float> gp0(nzpad * nxpad, 0);
  std::vector<float> gp1(nzpad * nxpad, 0);


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
}

void calMaxAlpha2_3(const Velocity &exvel,  const float *grad, float dt, float dx, float maxdv,
                    float &ret_alpha2, float &ret_alpha3) {
  const int nx = exvel.nx;
  const int nz = exvel.nz;

  const std::vector<float> &vel = exvel.dat;
  float alpha2 = FLT_MAX;
  for (int i = 0; i < nx * nz; i++) {
    float tmpv = dx / (dt * std::sqrt(vel[i]));
    tmpv -= maxdv;
    tmpv = (dx / (dt * tmpv)) * (dx / (dt * tmpv));
    if (std::fabs(grad[i]) < 1e-10 ) {
      continue;
    }
    if (alpha2 > (tmpv - vel[i]) / std::fabs(grad[i])) {
      alpha2 = (tmpv - vel[i]) / std::fabs(grad[i]);
    }
  }

  /// return the value
  ret_alpha2 = alpha2;
  ret_alpha3 =  2 * alpha2;
}

void update_vel(const float *vel, const float *grad, float size, float steplen, float vmin, float vmax, float *new_vel) {
  if (vmax <= vmin) {
    ERROR() << format("vmax(%f) < vmin(%f)") % vmax % vmin;
    exit(0);
  }

  for (int i = 0; i < size; i++) {
    new_vel[i] = vel[i] + steplen * grad[i];
    if (new_vel[i] > vmax) {
      new_vel[i] = vmax;
    }
    if (new_vel[i] < vmin) {
      new_vel[i] = vmin;
    }
  }
}

void initAlpha2_3(int ivel, float max_alpha3, float &initAlpha2, float &initAlpha3) {
  const float minAlpha   = 1.0E-7;
  const float resetAlpha = 1.0E-4;

  if (!PreservedAlpha::instance().getIsInit()[ivel]) {
    PreservedAlpha::instance().getIsInit()[ivel] = true;
    PreservedAlpha::instance().getAlpha()[ivel] = max_alpha3;
  }

  initAlpha3 = PreservedAlpha::instance().getAlpha()[ivel];
  initAlpha3 = initAlpha3 < minAlpha ? resetAlpha : initAlpha3;
  initAlpha2 = initAlpha3 * 0.5;
}

int calculate_obj_val(const Damp4t10d &fmMethod,
    const std::vector<float> &encsrc, const std::vector<float> &encobs,
    const float *grad,
    float vmin, float vmax, float steplen, float *obj_val_out) {

  int nx = fmMethod.getVelocity().nx;
  int nz = fmMethod.getVelocity().nz;
  int nt = fmMethod.getnt();

  int size = nx *  nz;

  const float *vel = &fmMethod.getVelocity().dat[0];
  float *new_vel = (float *)malloc(sizeof(float) * size);
  update_vel(vel, grad, size, steplen, vmin, vmax, new_vel);

  Damp4t10d updateMethod = fmMethod;
  Velocity updateVel(std::vector<float>(new_vel, new_vel + size), nx, nz);
  updateMethod.bindVelocity(updateVel);


  //forward modeling
  int ng = fmMethod.getng();
  std::vector<float> dcal(nt * ng);
  forwardModeling(updateMethod, encsrc, dcal, nt);


  updateMethod.removeDirectArrival(&dcal[0]);

  std::vector<float> vdiff(nt * ng, 0);
  vectorMinus(encobs, dcal, vdiff);
  float val = cal_objective(&vdiff[0], vdiff.size());

  DEBUG() << format("curr_alpha = %e, pure object value = %e") % steplen % val;

  *obj_val_out = val;

  return 0;
}

void selectAlpha(const Damp4t10d &fmMethod,
    const std::vector<float> &encsrc, const std::vector<float> &encobs, const float *grad,
    float obj_val1, float vmin, float vmax, float maxAlpha3,
    float &_alpha2, float &_obj_val2, float &_alpha3, float &_obj_val3, bool &toParabolicFit) {
  TRACE() << "SELECTING THE RIGHT OBJECTIVE VALUE 3";

  float alpha3 = _alpha3;
  float alpha2 = _alpha2;
  float obj_val2, obj_val3 = 0;

  calculate_obj_val(fmMethod, encsrc, encobs, grad, vmin, vmax, alpha2, &obj_val2);
  calculate_obj_val(fmMethod, encsrc, encobs, grad, vmin, vmax, alpha3, &obj_val3);

  DEBUG() << "BEFORE TUNNING";
  DEBUG() << __FUNCTION__ << format(" alpha1 = %e, obj_val1 = %e") % 0. % obj_val1;
  DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % alpha2 % obj_val2;
  DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % alpha3 % obj_val3;


  TRACE() << "maintain a set to store alpha2 that we ever tuned";
  std::set<ParaPoint, bool (*)(const ParaPoint &, const ParaPoint &) > tunedAlpha(parabolicLessComp);
  tunedAlpha.insert(std::make_pair(alpha2, obj_val2));

  DEBUG() << "BEGIN TUNING";
  /// obj_val2 might be quite large, so we should make it smaller by halfing alpha2

  int iter = 0;
  for (; iter < max_iter_select_alpha3 && obj_val2 > obj_val1; iter++) {

    /// pass the property of alpha2 to alpha3
    alpha3 = alpha2;
    obj_val3 = obj_val2;

    /// update alpha2
    alpha2 /= 2;
    calculate_obj_val(fmMethod, encsrc, encobs, grad, vmin, vmax, alpha2, &obj_val2);

    /// store it
    tunedAlpha.insert(std::make_pair(alpha2, obj_val2));
    DEBUG() << __FUNCTION__ << format(" iter = %d, alpha2 = %e, obj_val2 = %e") % iter % alpha2 % obj_val2;
    DEBUG() << __FUNCTION__ << format(" iter = %d, alpha3 = %e, obj_val3 = %e\n") % iter % alpha3 % obj_val3;
  }

  DEBUG() << "SELECT A BETTER ALPHA2 IN " << iter << " ITERS";


  DEBUG() << "tunedAlpha size: " << tunedAlpha.size();
  for (std::set<ParaPoint, bool (*)(const ParaPoint &, const ParaPoint &) >::iterator it = tunedAlpha.begin();
      it != tunedAlpha.end(); ++it) {
    DEBUG() << format("alpha %e, obj %e") % it->first % it->second;
  }

  TRACE() << "check if we need to forward tuning";
  TRACE() << "after halfing in the previous step, obj_val2 might still be larger than obj_val1"
          "then we should stop tunting and choose a best alpha2 ever got";
  if (obj_val2 > obj_val1) {
    DEBUG() << "UNABLE TO TUNING A ALPHA2 BY HALFING";
    DEBUG() << "SELECT A BEST ALPHA2 EVER GOT";
    std::set<ParaPoint>::iterator it = tunedAlpha.begin();
    _alpha2 = it->first;
    _obj_val2 = it->second;

    _alpha3 = std::min(_alpha2 * 2, maxAlpha3);
    calculate_obj_val(fmMethod, encsrc, encobs, grad, vmin, vmax, alpha3, &_obj_val3);

    toParabolicFit = false;

    DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
    DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;
    return;
  }

  TRACE() << "now we can make sure that obj_val2 < obj_val1";

  const float alpha1 = 0;
  float linearFitAlph3 = (obj_val2 - obj_val1) / (alpha2 - alpha1) * (alpha3 - alpha1) + obj_val1;
  DEBUG() << __FUNCTION__ << format(" linear fit alpha3 = %e ") % linearFitAlph3;

  TRACE() << "keep the alpha we tuned";
  tunedAlpha.clear();
  tunedAlpha.insert(std::make_pair(alpha3, obj_val3));

  while (obj_val3 < linearFitAlph3 && obj_val3 < obj_val1 && alpha3 < maxAlpha3) {
    TRACE() << "if in this case, we should enlarge alpha3";
    alpha2 = alpha3;
    obj_val2 = obj_val3;

    alpha3 = std::min(alpha3 * 2, maxAlpha3);
    calculate_obj_val(fmMethod,  encsrc, encobs, grad, vmin, vmax, alpha3, &obj_val3);

    tunedAlpha.insert(std::make_pair(alpha3, obj_val3));

    DEBUG() << __FUNCTION__ << format(" tune alpha3, alpha2 = %e, obj_val2 = %e") % alpha2 % obj_val2;
    DEBUG() << __FUNCTION__ << format(" tune alpha3, alpha3 = %e, obj_val3 = %e") % alpha3 % obj_val3;
  }


  TRACE() << "If we couldnot tune a good alpha3";
  if (alpha3 > maxAlpha3 + 0.1) {
    DEBUG() << "UNABLE TO TUNING A ALPHA3 BY DOUBLING";
    DEBUG() << "SELECT A BEST ALPHA3 EVER GOT";
    std::set<ParaPoint>::iterator it = tunedAlpha.begin();
    _alpha3 = it->first;
    _obj_val3 = it->second;

    _alpha2 = _alpha3 / 2;
    calculate_obj_val(fmMethod,  encsrc, encobs, grad, vmin, vmax, alpha2, &_obj_val2);

    toParabolicFit = false;

    DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
    DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;
    return;
  }

  /// return objval2 and objval3
  toParabolicFit = true;
  _alpha2 = alpha2;
  _alpha3 = alpha3;
  _obj_val2 = obj_val2;
  _obj_val3 = obj_val3;

  DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
  DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;
}




float calStepLen(const Damp4t10d &fmMethod,
    const std::vector<float> &encsrc, const std::vector<float> &encobs,
    const std::vector<float> &updateDirection, int iter, int ivel,
    float obj_val1, float min_vel, float max_vel) {

  float dt = fmMethod.getdt();
  float dx = fmMethod.getdx();

  TRACE() << "Calcuate step length";
  TRACE() << "calculate the initial value of alpha2 and alpha3";
  float max_alpha2, max_alpha3;
  calMaxAlpha2_3(fmMethod.getVelocity(), &updateDirection[0], dt, dx, maxdv, max_alpha2, max_alpha3);
  DEBUG() << format("               max_alpha2 = %e,  max_alpha3: = %e") % max_alpha2 % max_alpha3;

  float alpha1 = 0, alpha2, alpha3;
  initAlpha2_3(ivel, max_alpha3, alpha2, alpha3);
  DEBUG() << format("after init alpha,  alpha2 = %e,      alpha3: = %e") % alpha2 % alpha3;

  float obj_val2, obj_val3;
  bool toParabolic;

  selectAlpha(fmMethod, encsrc, encobs, &updateDirection[0], obj_val1, min_vel, max_vel, max_alpha3, alpha2, obj_val2, alpha3, obj_val3, toParabolic);

  float alpha4, obj_val4;
  if (toParabolic) {
    DEBUG() << "parabolic fit";
    parabolaVertex(alpha1, obj_val1, alpha2, obj_val2, alpha3, obj_val3, max_alpha3, alpha4, obj_val4);
    if (alpha4 > max_alpha3) {
      DEBUG() << format("alpha4 = %e, max_alpha3 = %e") % alpha4 % max_alpha3;
      DEBUG() << format("alpha4 is greater than max_alpha3, set it to alpha3");
      alpha4 = max_alpha3;
    }
  } else {
    DEBUG() << "NO need to perform parabolic fit";
    alpha4 = alpha3;
    obj_val4 = obj_val3;
  }

  INFO() << format("In calculate_steplen(): iter %d  alpha  = %e total obj_val1 = %e") % iter % alpha1 % obj_val1;
  INFO() << format("In calculate_steplen(): iter %d  alpha2 = %e total obj_val2 = %e") % iter % alpha2 % obj_val2;
  INFO() << format("In calculate_steplen(): iter %d  alpha3 = %e total obj_val3 = %e") % iter % alpha3 % obj_val3;
  INFO() << format("In calculate_steplen(): iter %d  alpha4 = %e total obj_val4 = %e\n") % iter % alpha4 % obj_val4;

  PreservedAlpha::instance().getAlpha()[ivel] = alpha4;
  return alpha4;
}

}

EssFwiFramework::EssFwiFramework(Damp4t10d &method, const std::vector<float> &_wlt,
    const std::vector<float> &_dobs) :
    fmMethod(method), wlt(_wlt), dobs(_dobs),
    ns(method.getns()), ng(method.getng()), nt(method.getnt()),
    nx(method.getnx()), nz(method.getnz()), dx(method.getdx()), dt(method.getdt())
//    ,g0(nx*nz, 0), updateDirection(nx*nz, 0)
{
  g0.resize(nx*nz, 0);
  updateDirection.resize(nx*nz, 0);
}

void EssFwiFramework::epoch(int iter, int ivel) {

  // create random codes
  const std::vector<int> encodes = RandomCode::genPlus1Minus1(ns);
  std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

  Encoder encoder(encodes);
  std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
  std::vector<float> encsrc  = encoder.encodeSource(wlt);

  {
    char buf[BUFSIZ];
    sprintf(buf, "encobs%d.rsf", iter);
//    sfFloatWrite2d(buf, &encobs[0], nt, ng);

    sprintf(buf, "encsrc%d.rsf", iter);
//    sfFloatWrite1d(buf, &encsrc[0], encsrc.size());

    DEBUG() << format("sum wlt %.20f") % sum(wlt);
//    sprintf(buf, "exvel%d.rsf", iter);
    //sfFloatWrite2d(buf, &exvel.dat[0], exvel.nz, exvel.nx);
  }
//
  std::vector<float> dcal(nt * ng, 0);
  forwardModeling(fmMethod, encsrc, dcal, nt);
//
  {
    char buf[BUFSIZ];
    sprintf(buf, "calobs%d.rsf", iter);
    //sfFloatWrite2d(buf, &dcal[0], ng, nt);
  }
//
  fmMethod.removeDirectArrival(&encobs[0]);
  fmMethod.removeDirectArrival(&dcal[0]);
//
  {
    char buf[BUFSIZ];
    sprintf(buf, "rmdcalobs%d.rsf", iter);
    //sfFloatWrite2d(buf, &dcal[0], ng, nt);
  }
//
  std::vector<float> vsrc(nt * ng, 0);
  vectorMinus(encobs, dcal, vsrc);
  float obj1 = cal_objective(&vsrc[0], vsrc.size());
  DEBUG() << format("obj: %e") % obj1;
//  exit(0);
//
  transVsrc(vsrc, nt, ng, dt);
//
//
  std::vector<float> g1(nx * nz, 0);
  hello2(fmMethod, encsrc, vsrc, g1, nt, dt);
//  {
//    char buf[BUFSIZ];
//    sprintf(buf, "grad%d.rsf", iter);
//    //sfFloatWrite2d(buf, &g1[0], exvel.nz, exvel.nx);
//  }
////    exit(0);
//
  fmMethod.maskGradient(&g1[0]);
//  {
//    char buf[BUFSIZ];
//    sprintf(buf, "mgrad%d.rsf", iter);
//    //sfFloatWrite2d(buf, &g1[0], exvel.nz, exvel.nx);
//  }
//
//
//  {
//    char buf[BUFSIZ];
//    sprintf(buf, "pre%d.rsf", iter);
//    //sfFloatWrite2d(buf, &g0[0], exvel.nz, exvel.nx);
//  }

//  prevCurrCorrDirection(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);
    gradUpdator(&g0[0], &g1[0], &updateDirection[0], g0.size()); /// for enfwi
//
//  {
//    char buf[BUFSIZ];
//    sprintf(buf, "g0%d.rsf", iter);
//    //sfFloatWrite2d(buf, &g0[0], exvel.nz, exvel.nx);
//
//    sprintf(buf, "update%d.rsf", iter);
//    //sfFloatWrite2d(buf, &updateDirection[0], exvel.nz, exvel.nx);
//  }
//
  float min_vel = (dx / dt / vmax) * (dx / dt / vmax);
  float max_vel = (dx / dt / vmin) * (dx / dt / vmin);




  DEBUG() << format("vmax: %f, vmin: %f, minv: %f, maxv: %f") % vmax % vmin % min_vel % max_vel;
  float steplen = calStepLen(fmMethod, encsrc, encobs, updateDirection, iter, ivel, obj1, min_vel, max_vel);






  Velocity &exvel = fmMethod.getVelocity();
  TRACE() << "Update velocity model";
  update_vel(&exvel.dat[0], &updateDirection[0], exvel.dat.size(), steplen, min_vel, max_vel, &exvel.dat[0]);
////    sfFloatWrite2d("updatevel.rsf", &exvel.dat[0], exvel.nz, exvel.nx);
//
  fmMethod.refillBoundary(&exvel.dat[0]);
////    sfFloatWrite2d("updatevel-refilled.rsf", &exvel.dat[0], exvel.nz, exvel.nx);
}

void EssFwiFramework::writeVel(sf_file file) const {
  fmMethod.sfWriteVel(file);
}
