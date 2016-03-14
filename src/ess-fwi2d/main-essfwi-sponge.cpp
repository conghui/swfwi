
/* Time domain full waveform inversion
Note: This serial FWI is merely designed to help the understanding of
beginners. Enquist absorbing boundary condition (A2) is applied!
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Important references:
    [1] Clayton, Robert, and Björn Engquist. "Absorbing boundary
  conditions for acoustic and elastic wave equations." Bulletin
  of the Seismological Society of America 67.6 (1977): 1529-1540.
    [2] Tarantola, Albert. "Inversion of seismic reflection data in the
  acoustic approximation." Geophysics 49.8 (1984): 1259-1266.
    [3] Pica, A., J. P. Diet, and A. Tarantola. "Nonlinear inversion
  of seismic reflection data in a laterally invariant medium."
  Geophysics 55.3 (1990): 284-292.
    [4] Dussaud, E., Symes, W. W., Williamson, P., Lemaistre, L.,
  Singer, P., Denel, B., & Cherrett, A. (2008). Computational
  strategies for reverse-time migration. In SEG Technical Program
  Expanded Abstracts 2008 (pp. 2267-2271).
    [5] Hager, William W., and Hongchao Zhang. "A survey of nonlinear
  conjugate gradient methods." Pacific journal of Optimization
  2.1 (2006): 35-58.
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
#include "zjh4t10dsponge.h"
#include "sfutil.h"
#include "aux.h"
#include "preserved-alpha.h"

static const int max_iter_select_alpha3 = 5;
static const float vmax = 5500;
static const float vmin = 1500;
static const float maxdv = 200;
typedef std::pair<float, float> ParaPoint;

bool parabolicLessComp(const ParaPoint &a, const ParaPoint &b) {
  return a.second - b.second < 1e-10;
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

    DEBUG() << format("beta %e") % beta;
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

void forwardModeling(const Zjh4t10dSponge &fmMethod,
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

void hello2(const Zjh4t10dSponge &fmMethod,
    const ShotPosition &allSrcPos, const std::vector<float> &encSrc,
    const ShotPosition &allGeoPos, const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt)
{
  int nxpad = fmMethod.getVelocity().nx;
  int nzpad = fmMethod.getVelocity().nz;
  int ns = allSrcPos.ns;
  int ng = allGeoPos.ns;

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

int calculate_obj_val(const Zjh4t10dSponge &fmMethod,
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

  Zjh4t10dSponge updateMethod = fmMethod;
  Velocity updateVel(std::vector<float>(new_vel, new_vel + size), nx, nz);
  updateMethod.bindVelocity(updateVel);

//  sfFloatWrite2d("upvel.rsf", &updateMethod.getVelocity().dat[0], nz, nx);

  //forward modeling
  int ng = fmMethod.getng();
  std::vector<float> dcal(nt * ng);
  forwardModeling(updateMethod, encsrc, dcal, nt);

//  sfFloatWrite2d("11dcal.rsf", &dcal[0], ng, nt);

  updateMethod.removeDirectArrival(&dcal[0]);

//  sfFloatWrite2d("22dcal.rsf", &dcal[0], ng, nt);

  std::vector<float> vdiff(nt * ng, 0);
  vectorMinus(encobs, dcal, vdiff);
  float val = cal_objective(&vdiff[0], vdiff.size());

  DEBUG() << format("curr_alpha = %e, pure object value = %e") % steplen % val;

  *obj_val_out = val;

  return 0;
}

void selectAlpha(const Zjh4t10dSponge &fmMethod,
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

//  exit(0);

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
//    calculate_obj_val(dim, config, shot, grad, vel, vmin, vmax, alpha2, &obj_val2);

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
//    calculate_obj_val(dim, config, shot, grad, vel, vmin, vmax, _alpha3, &_obj_val3);

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
//    calculate_obj_val(dim, config, shot, grad, vel, vmin, vmax, alpha3, &obj_val3);

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
//    calculate_obj_val(dim, config, shot, grad, vel, vmin, vmax, _alpha2, &_obj_val2);

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

void calcParabolaVertex(float x1, float y1, float x2, float y2, float x3, float y3, float &xv, float &yv) {
  double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
  double A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
  double B     = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
  double C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

  xv = -B / (2 * A);
  yv = C - B * B / (4 * A);
}

void calcParabolaVertexEnhanced(float x1, float y1, float x2, float y2, float x3, float y3, float max_alpha3, float &xv, float &yv) {
  double k2 = (y3 - y2) / (x3 - x2);
  double k1 = (y2 - y1) / (x2 - x1);

  calcParabolaVertex(x1, y1, x2, y2, x3, y3, xv, yv);

  if (std::abs(k2 - k1) < 0.001 * (std::max(std::abs(k2), std::abs(k1))) ||
      (xv == -std::numeric_limits<double>::quiet_NaN())) {
    WARNING() << "THE SET OF POINTS DON'T FIT PARABOLIC WELL, SET y TO -NAN ON PURPOSE JUST FOR INDICATION";
    xv = std::min(2 * x3, max_alpha3);
    yv = -std::numeric_limits<double>::quiet_NaN(); /// indicating what's happening
  }
}


float calStepLen(const Zjh4t10dSponge &fmMethod,
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
    calcParabolaVertexEnhanced(alpha1, obj_val1, alpha2, obj_val2, alpha3, obj_val3, max_alpha3, alpha4, obj_val4);
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


int main(int argc, char *argv[]) {

  /* initialize Madagascar */
  sf_init(argc, argv);

  Logger::instance().init("essfwi");

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

  Zjh4t10dSponge fmMethod(allSrcPos, allGeoPos, dt, dx, fm, nb, nt);

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);
  Velocity exvel = fmMethod.expandDomain(v0);

  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);


  std::vector<float> dobs(ns * nt * ng);     /* all observed data */
  std::vector<float> g0(exvel.nx * exvel.nz, 0); /* gradient at previous step */

  ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);
  //sfFloatWrite1d("orgdata.rsf", &dobs[0], ns * nt * ng);

  std::vector<float> updateDirection(exvel.nx * exvel.nz, 0);

  for (int iter = 0; iter < params.niter; iter++) {
    boost::timer::cpu_timer timer;

    // create random codes
    const std::vector<int> encodes = RandomCode::genPlus1Minus1(params.ns);
//    std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

    Encoder encoder(encodes);
    std::vector<float> encobs = encoder.encodeObsData(dobs, params.nt, params.ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    {
      char buf[BUFSIZ];
      sprintf(buf, "encobs%d.rsf", iter);
      //sfFloatWrite2d(buf, &encobs[0], nt, ng);

      sprintf(buf, "encsrc%d.rsf", iter);
      //sfFloatWrite1d(buf, &encsrc[0], encsrc.size());

      sprintf(buf, "exvel%d.rsf", iter);
      //sfFloatWrite2d(buf, &exvel.dat[0], exvel.nz, exvel.nx);
    }

    std::vector<float> dcal(nt * ng, 0);
    forwardModeling(fmMethod, encsrc, dcal, nt);

    {
      char buf[BUFSIZ];
      sprintf(buf, "calobs%d.rsf", iter);
      //sfFloatWrite2d(buf, &dcal[0], ng, nt);
    }

    fmMethod.removeDirectArrival(&encobs[0]);
    fmMethod.removeDirectArrival(&dcal[0]);

    {
      char buf[BUFSIZ];
      sprintf(buf, "rmdcalobs%d.rsf", iter);
      //sfFloatWrite2d(buf, &dcal[0], ng, nt);
    }

    std::vector<float> vsrc(nt * ng, 0);
    vectorMinus(encobs, dcal, vsrc);
    float obj1 = cal_objective(&vsrc[0], vsrc.size());
    DEBUG() << format("obj: %e") % obj1;
//    exit(0);

    transVsrc(vsrc, nt, ng, dt);

//    forwardPropagate(fmMethod, allSrcPos, encsrc, nt);

    std::vector<float> g1(exvel.nx * exvel.nz, 0);
    hello2(fmMethod, allSrcPos, encsrc, allGeoPos, vsrc, g1, nt, dt);
    {
      char buf[BUFSIZ];
      sprintf(buf, "grad%d.rsf", iter);
      //sfFloatWrite2d(buf, &g1[0], exvel.nz, exvel.nx);
    }
//    exit(0);

    fmMethod.maskGradient(&g1[0]);
    {
      char buf[BUFSIZ];
      sprintf(buf, "mgrad%d.rsf", iter);
      //sfFloatWrite2d(buf, &g1[0], exvel.nz, exvel.nx);
    }


    {
      char buf[BUFSIZ];
      sprintf(buf, "pre%d.rsf", iter);
      //sfFloatWrite2d(buf, &g0[0], exvel.nz, exvel.nx);
    }
    prevCurrCorrDirection(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);

    {
      char buf[BUFSIZ];
      sprintf(buf, "g0%d.rsf", iter);
      //sfFloatWrite2d(buf, &g0[0], exvel.nz, exvel.nx);

      sprintf(buf, "update%d.rsf", iter);
      //sfFloatWrite2d(buf, &updateDirection[0], exvel.nz, exvel.nx);
    }

    const int ivel = 0;
    float min_vel = (dx / dt / vmax) * (dx / dt / vmax);
    float max_vel = (dx / dt / vmin) * (dx / dt / vmin);

    DEBUG() << format("vmax: %f, vmin: %f, minv: %f, maxv: %f") % vmax % vmin % min_vel % max_vel;
    float steplen = calStepLen(fmMethod, encsrc, encobs, updateDirection, iter, ivel, obj1, min_vel, max_vel);

    TRACE() << "Update velocity model";
    update_vel(&exvel.dat[0], &updateDirection[0], exvel.dat.size(), steplen, min_vel, max_vel, &exvel.dat[0]);
//    sfFloatWrite2d("updatevel.rsf", &exvel.dat[0], exvel.nz, exvel.nx);

    fmMethod.refillBoundary(&exvel.dat[0]);
//    sfFloatWrite2d("updatevel-refilled.rsf", &exvel.dat[0], exvel.nz, exvel.nx);
    fmMethod.sfWriteVel(params.vupdates);
  } /// end of iteration

  sf_close();

  return 0;
}
