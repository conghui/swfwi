
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
#include "damp4t10d.h"
#include "sfutil.h"
#include "aux.h"


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

float cal_obj_derr_illum_grad(const EssFwiParams &params,
    float *derr,  /* output */
    float *illum, /* output */
    float *g1,    /* output */
    const float *vv,
    const float *encsrc,
    const float *encdobs,
    const int *sxz,
    const int *gxz)
{
  int nt = params.nt;
  int nz = params.nz;
  int nx = params.nx;
  int ng = params.ng;
  int ns = params.ns;
  float dt = params.dt;
  float dx = params.dx;
  float dz = params.dz;

  std::vector<float> bndr(nt * (2 * nz + nx), 0); /* boundaries for wavefield reconstruction */
  std::vector<float> dcal(ng, 0); /* calculated/synthetic seismic data */

  std::vector<float> sp0(nz * nx); /* source wavefield p0 */
  std::vector<float> sp1(nz * nx); /* source wavefield p1 */
  std::vector<float> sp2(nz * nx); /* source wavefield p2 */
  std::vector<float> gp0(nz * nx); /* geophone/receiver wavefield p0 */
  std::vector<float> gp1(nz * nx); /* geophone/receiver wavefield p1 */
  std::vector<float> gp2(nz * nx); /* geophone/receiver wavefield p2 */
  std::vector<float> lap(nz * nx); /* laplace of the source wavefield */

  float dtx = dt / dx;
  float dtz = dt / dz;
  int nb = 0; // no expanded boundary

  std::fill(sp0.begin(), sp0.end(), 0);
  std::fill(sp1.begin(), sp1.end(), 0);

  for (int it = 0; it < nt; it++) {
    // add all sources to the wave field
    add_source(&sp1[0], &encsrc[it * ns], &sxz[0], ns, nz, nb, true);
    step_forward(&sp0[0], &sp1[0], &sp2[0], vv, dtz, dtx, nz, nx);
    cycleSwap(sp0, sp1, sp2); // cycle swap
    rw_bndr(&bndr[it * (2 * nz + nx)], &sp0[0], nz, nx, true);
    record_seis(&dcal[0], gxz, &sp0[0], ng, nz, nb);
    cal_residuals(&dcal[0], &encdobs[it * ng], &derr[it * ng], ng);
  }

  std::swap(sp0, sp1);

  std::fill(gp0.begin(), gp0.end(), 0);
  std::fill(gp1.begin(), gp1.end(), 0);

  for (int it = nt - 1; it > -1; it--) {
    rw_bndr(&bndr[it * (2 * nz + nx)], &sp1[0], nz, nx, false);
    step_backward(illum, &lap[0], &sp0[0], &sp1[0], &sp2[0], vv, dtz, dtx, nz, nx);
    // remove all source from wave field
    add_source(&sp1[0], &encsrc[it * ns], &sxz[0], ns, nz, nb, false);

    add_source(&gp1[0], &derr[it * ng], gxz, ng, nz, nb, true);

    step_forward(&gp0[0], &gp1[0], &gp2[0], vv, dtz, dtx, nz, nx);

    cal_gradient(&g1[0], &lap[0], &gp1[0], nz, nx);

    cycleSwap(sp0, sp1, sp2);
    cycleSwap(gp0, gp1, gp2);
  }

  float obj = cal_objective(&derr[0], ng * nt);

  return obj;
}

float calVelUpdateStepLen(const EssFwiParams &params,
    const float *vtmp,
    const float *encsrc,
    const float *encdobs,
    const int *sxz,
    const int *gxz,
    const float *derr,
    float epsil
    )
{
  int nt = params.nt;
  int nz = params.nz;
  int nx = params.nx;
  int ng = params.ng;
  int ns = params.ns;
  float dt = params.dt;
  float dx = params.dx;
  float dz = params.dz;

  std::vector<float> dcal(ng, 0); /* calculated/synthetic seismic data */
  std::vector<float> sp0(nz * nx); /* source wavefield p0 */
  std::vector<float> sp1(nz * nx); /* source wavefield p1 */
  std::vector<float> sp2(nz * nx); /* source wavefield p2 */

  std::vector<float> alpha1(ng, 0); /* numerator of alpha, length=ng */
  std::vector<float> alpha2(ng, 0); /* denominator of alpha, length=ng */

  float dtx = dt / dx;
  float dtz = dt / dz;
  int nb = 0; // no expanded boundary

  std::fill(sp0.begin(), sp0.end(), 0);
  std::fill(sp1.begin(), sp1.end(), 0);
  for (int it = 0; it < nt; it++) {
    // add all sources to the wave field
    add_source(&sp1[0], &encsrc[it * ns], &sxz[0], ns, nz, nb, true);
    step_forward(&sp0[0], &sp1[0], &sp2[0], &vtmp[0], dtz, dtx, nz, nx);

    std::swap(sp0, sp1);
    std::swap(sp1, sp2);

    record_seis(&dcal[0], gxz, &sp0[0], ng, nz, nb);
    sum_alpha12(&alpha1[0], &alpha2[0], &dcal[0], &encdobs[it * ng], &derr[it * ng], ng);
  }

  float alpha = cal_alpha(&alpha1[0], &alpha2[0], epsil, ng);

  return alpha;
}

void forwardModeling(const Damp4t10d &fmMethod,
    const ShotPosition &allSrcPos, const ShotPosition &allGeoPos,
    const std::vector<float> &encSrc, std::vector<float> &dobs,
    int nt)
{
    int nxpad = fmMethod.getVelocity().nx;
    int nzpad = fmMethod.getVelocity().nz;
    int ns = allSrcPos.ns;
    int ng = allGeoPos.ns;

    boost::timer::cpu_timer timer;

    std::vector<float> p0(nzpad * nxpad, 0);
    std::vector<float> p1(nzpad * nxpad, 0);

    for(int it=0; it<nt; it++) {

      fmMethod.addSource(&p1[0], &encSrc[it * ns], allSrcPos);

      fmMethod.stepForward(&p0[0], &p1[0]);

      fmMethod.recordSeis(&dobs[it*ng], &p0[0], allGeoPos);

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

  sfFloatWrite2d("vsrc.rsf", &trans[0], nt, ng);

  matrix_transpose(&trans[0], &vsrc[0], nt, ng);
}

void forwardPropagate(const Damp4t10d &fmMethod,
    const ShotPosition &allSrcPos, const std::vector<float> &encSrc,
    int nt)
{
  const int check_step = 5;

  int nxpad = fmMethod.getVelocity().nx;
  int nzpad = fmMethod.getVelocity().nz;
  int ns = allSrcPos.ns;

  std::vector<float> p0(nzpad * nxpad, 0);
  std::vector<float> p1(nzpad * nxpad, 0);

  for(int it=0; it<nt; it++) {
    fmMethod.addSource(&p1[0], &encSrc[it * ns], allSrcPos);
    fmMethod.stepForward(&p0[0], &p1[0]);
    std::swap(p1, p0);

    if ((it > 0) && (it != (nt - 1)) && !(it % check_step)) {
      char check_file_name1[64];
      char check_file_name2[64];
      const char *checkPointDir = std::getenv("CHECKPOINTDIR");
      sprintf(check_file_name1, "%s/check_time_%d_1.su", checkPointDir, it);
      sprintf(check_file_name2, "%s/check_time_%d_2.su", checkPointDir, it);
      writeBin(std::string(check_file_name1), &p0[0], p0.size() * sizeof(float));
      writeBin(std::string(check_file_name2), &p1[0], p1.size() * sizeof(float));
    }
  }

  char check_file_name1[64];
  char check_file_name2[64];
  const char *checkPointDir = std::getenv("CHECKPOINTDIR");
  sprintf(check_file_name1, "%s/check_time_last_1.su", checkPointDir);
  sprintf(check_file_name2, "%s/check_time_last_2.su", checkPointDir);
  writeBin(std::string(check_file_name1), &p0[0], p0.size() * sizeof(float));
  writeBin(std::string(check_file_name2), &p1[0], p1.size() * sizeof(float));
}

static void cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale) {
  int i = 0;
  for (i = 0; i < model_size; i ++) {
    image[i] -= src_wave[i] * vsrc_wave[i] * scale;
  }

}


void hello(const Damp4t10d &fmMethod,
    const ShotPosition &allSrcPos, const std::vector<float> &encSrc,
    const ShotPosition &allGeoPos, const std::vector<float> &vsrc,
    std::vector<float> &g0,
    int nt, float dt)
{
  const int check_step = 50;

  int nxpad = fmMethod.getVelocity().nx;
  int nzpad = fmMethod.getVelocity().nz;
  int ns = allSrcPos.ns;
  int ng = allGeoPos.ns;

  std::vector<float> sp0(nzpad * nxpad, 0);
  std::vector<float> sp1(nzpad * nxpad, 0);
  std::vector<float> gp0(nzpad * nxpad, 0);
  std::vector<float> gp1(nzpad * nxpad, 0);

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
      printf("reading %s and %s\n", check_file_name1, check_file_name2);
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

    fmMethod.subSource(&sp0[0], &encSrc[it * ns], allSrcPos);
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
//    }

    fmMethod.stepBackward(&gp0[0], &gp1[0]);
//    {
//      char buf[256];
//      sprintf(buf, "gp0aftfm%d.rsf", it);
//      sfFloatWrite2d(buf, &gp0[0], nzpad, nxpad);
//    }

    std::swap(gp1, gp0);

//    if (it == 0) {
//      char buf[256];
//      sprintf(buf, "sfield%d.rsf", it);
//      sfFloatWrite2d(buf, &sp0[0], nzpad, nxpad);
//
//      sprintf(buf, "vfield%d.rsf", it);
//      sfFloatWrite2d(buf, &gp0[0], nzpad, nxpad);
////    }
//
//      if (it == 397) exit(0);

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

void calStepLen(float dt, float dx) {
  TRACE() << "Calcuate step length";
  const float vmax = 5500;
  const float vmin = 1500;
  float min_vel = (dx / dt / vmax) * (dx / dt / vmax);
  float max_vel = (dx / dt / vmin) * (dx / dt / vmin);
  DEBUG() << format("vmax: %f, vmin: %f, minv: %f, maxv: %f") % vmax % vmin % min_vel % max_vel;

//  DEBUG() << format
  TRACE() << "calculate the initial value of alpha2 and alpha3";
  float max_alpha2, max_alpha3;
//  calMaxAlpha2_3(dim, config, grad, vel, max_alpha2, max_alpha3)

//  float min_vel = (dim.dx / config.dt / config.velocity_max) * (dim.dx / config.dt / config.velocity_max);
//  float max_vel = (dim.dx / config.dt / config.velocity_min) * (dim.dx / config.dt / config.velocity_min);
}
//float subSquareSum(const std::vector<float> &a, const std::vector<float> &b) {
//  float sum = 0;
//  for (size_t i = 0; i < a.size(); i++) {
//    sum += (a[i] - b[i]) * (a[i] - b[i]);
//  }
//
//  return sum;
//}

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

  // set random seed
  const int seed = 10;
  srand(seed);

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);

  Damp4t10d fmMethod(dt, params.dx, nb);

  Velocity exvel = fmMethod.expandDomain(v0);
  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);

  std::vector<float> dobs(ns * nt * ng);     /* all observed data */
//  std::vector<float> cg(params.nz * params.nx, 0);    /* conjugate gradient */
  std::vector<float> g0(exvel.nx * exvel.nz, 0); /* gradient at previous step */
//  std::vector<float> objval(params.niter, 0); /* objective/misfit function */

  ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);
//  sfFloatWrite1d("orgdata.rsf", &dobs[0], ns * nt * ng);

  for (int iter = 0; iter < params.niter; iter++) {
    boost::timer::cpu_timer timer;
//    std::vector<float> g1(params.nz * params.nx, 0);    /* gradient at curret step */
//    std::vector<float> derr(params.ng * params.nt, 0); /* residual/error between synthetic and observation */
//    std::vector<float> illum(params.nz * params.nx, 0); /* illumination of the source wavefield */
//    std::vector<float> vtmp(params.nz * params.nx, 0);  /* temporary velocity computed with epsil */

    // create random codes
    const std::vector<int> encodes = RandomCode::genPlus1Minus1(params.ns);
    std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

    Encoder encoder(encodes);
    std::vector<float> encobs = encoder.encodeObsData(dobs, params.nt, params.ng);
    DEBUG() << format("sum encdobs %f") % sum(encobs);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    sfFloatWrite2d("encobs.rsf", &encobs[0], nt, ng);
    sfFloatWrite1d("encsrc.rsf", &encsrc[0], encsrc.size());

    std::vector<float> dcal(nt * ng, 0);
    forwardModeling(fmMethod, allSrcPos, allGeoPos, encsrc, dcal, nt);
    sfFloatWrite2d("calobs.rsf", &dcal[0], ng, nt);

    std::vector<float> vsrc(nt * ng, 0);
    vectorMinus(encobs, dcal, vsrc);
    float obj = cal_objective(&vsrc[0], vsrc.size());
    DEBUG() << format("obj: %e") % obj;

    transVsrc(vsrc, nt, ng, dt);

    forwardPropagate(fmMethod, allSrcPos, encsrc, nt);

    std::vector<float> g1(exvel.nx * exvel.nz, 0);
    hello(fmMethod, allSrcPos, encsrc, allGeoPos, vsrc, g1, nt, dt);
    sfFloatWrite2d("grad.rsf", &g1[0], exvel.nz, exvel.nx);

    fmMethod.maskGradient(&g1[0]);
    sfFloatWrite2d("mgrad.rsf", &g1[0], exvel.nz, exvel.nx);

    std::vector<float> updateDirection(exvel.nx * exvel.nz, 0);
    prevCurrCorrDirection(&g0[0], &g1[0], &updateDirection[0], g0.size(), iter);

    sfFloatWrite2d("g0.rsf", &g0[0], exvel.nz, exvel.nx);
    sfFloatWrite2d("update.rsf", &updateDirection[0], exvel.nz, exvel.nx);

    calStepLen(dt, params.dx);
    exit(0);
    /**
     * calculate local objective function & derr & illum & g1(gradient)
     */
//    float obj = cal_obj_derr_illum_grad(params, &derr[0], &illum[0], &g1[0], &vv[0], &encsrc[0], &encdobs[0], &sxz[0], &gxz[0]);

//    objval[iter] = iter == 0 ? obj0 = obj, 1.0 : obj / obj0;
//
//    float epsil = 0;
//    float beta = 0;
//    sf_floatwrite(&illum[0], params.nz * params.nx, params.illums);
//
//    scale_gradient(&g1[0], &vv[0], &illum[0], params.nz, params.nx, params.precon);
//    bell_smoothz(&g1[0], &illum[0], params.rbell, params.nz, params.nx);
//    bell_smoothx(&illum[0], &g1[0], params.rbell, params.nz, params.nx);
//    sf_floatwrite(&g1[0], params.nz * params.nx, params.grads);
//
//    beta = iter == 0 ? 0.0 : cal_beta(&g0[0], &g1[0], &cg[0], params.nz, params.nx);
//
//    cal_conjgrad(&g1[0], &cg[0], beta, params.nz, params.nx);
//    epsil = cal_epsilon(&vv[0], &cg[0], params.nz, params.nx);
//    cal_vtmp(&vtmp[0], &vv[0], &cg[0], epsil, params.nz, params.nx);
//
//    std::swap(g1, g0); // let g0 be the previous gradient
//
//    float alpha = calVelUpdateStepLen(params, &vtmp[0], &encsrc[0], &encdobs[0], &sxz[0], &gxz[0], &derr[0], epsil);
//
//    update_vel(&vv[0], &cg[0], alpha, params.nz, params.nx);
//
//    sf_floatwrite(&vv[0], params.nz * params.nx, params.vupdates);

    // output important information at each FWI iteration
//    INFO() << format("iteration %d obj=%f  beta=%f  epsil=%f  alpha=%f") % (iter + 1) % obj % beta % epsil % alpha;
//    INFO() << timer.format(2);

  } /// end of iteration

//  sf_floatwrite(&objval[0], params.niter, params.objs);

  sf_close();

  return 0;
}
