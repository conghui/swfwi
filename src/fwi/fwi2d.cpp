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

#include <mpi.h>
#include <time.h>
#include <algorithm>
#include <numeric>

#include <omp.h>

#include <boost/timer/timer.hpp>
#include "logger.h"
#include "global-params.h"

#include "common.h"

void cal_derr_illum_grad(const GlobalParams &params,
    float **vv,
    float *derr, /*output*/
    float **illum,
    float **g1,
    const float *wlt,
    const int *sxz,
    const int *gxz)
{
  std::vector<float> trans(params.ng * params.nt);
  std::vector<float> dobs(params.ng * params.nt);
  std::vector<float> bndr(params.nt * (2 * params.nz + params.nx), 0); /* boundaries for wavefield reconstruction */
  std::vector<float> dcal(params.ng, 0); /* calculated/synthetic seismic data */

  float **sp0 = sf_floatalloc2(params.nz, params.nx); /* source wavefield p0 */
  float **sp1 = sf_floatalloc2(params.nz, params.nx); /* source wavefield p1 */
  float **sp2 = sf_floatalloc2(params.nz, params.nx); /* source wavefield p2 */
  float **gp0 = sf_floatalloc2(params.nz, params.nx); /* geophone/receiver wavefield p0 */
  float **gp1 = sf_floatalloc2(params.nz, params.nx); /* geophone/receiver wavefield p1 */
  float **gp2 = sf_floatalloc2(params.nz, params.nx); /* geophone/receiver wavefield p2 */
  float **lap = sf_floatalloc2(params.nz, params.nx); /* laplace of the source wavefield */

  float dtx = params.dt / params.dx;
  float dtz = params.dt / params.dz;

  for (int is = 0; is < params.ns; is++) {
    sf_floatread(&trans[0], params.ng * params.nt, params.shots); /* Read local portion of input data */
    matrix_transpose(&trans[0], &dobs[0], params.nt, params.ng);

//      float sum_dobs = std::accumulate(dobs, dobs + params.nt * params.ng, 0.0f);
//      DEBUG() << format("iter %d, shot %d, sum_obs %.6f") % iter % is % sum_dobs;
    memset(sp0[0], 0, params.nz * params.nx * sizeof(float));
    memset(sp1[0], 0, params.nz * params.nx * sizeof(float));
    for (int it = 0; it < params.nt; it++) {
      add_source(sp1[0], &wlt[it], &sxz[is], 1, params.nz, true);
//      step_forward(sp0, sp1, sp2, vv, dtz, dtx, params.nz, params.nx);
      step_forward(sp0[0], sp1[0], sp2[0], vv[0], dtz, dtx, params.nz, params.nx);

      float **ptr = sp0;
      sp0 = sp1;
      sp1 = sp2;
      sp2 = ptr;

      rw_bndr(&bndr[it * (2 * params.nz + params.nx)], sp0[0], params.nz, params.nx, true);

      record_seis(&dcal[0], gxz, sp0[0], params.ng, params.nz);
      cal_residuals(&dcal[0], &dobs[it * params.ng], &derr[is * params.ng * params.nt + it * params.ng], params.ng);
    }

//    float sum_derr = std::accumulate(&derr[ik * params.ng * params.nt],
//        &derr[ik * params.ng * params.nt] + params.ng * params.nt, 0.0f);
//    DEBUG() << format("shot: %d, sum_derr %.6f") % (ik * size + rank) % sum_derr;
//
//
    float **ptr = sp0;
    sp0 = sp1;
    sp1 = ptr;

    memset(gp0[0], 0, params.nz * params.nx * sizeof(float));
    memset(gp1[0], 0, params.nz * params.nx * sizeof(float));

    for (int it = params.nt - 1; it > -1; it--) {
      rw_bndr(&bndr[it * (2 * params.nz + params.nx)], sp1[0], params.nz, params.nx, false);
      step_backward(illum[0], lap[0], sp0[0], sp1[0], sp2[0], vv[0], dtz, dtx, params.nz, params.nx);
      add_source(sp1[0], &wlt[it], &sxz[is], 1, params.nz, false);

      add_source(gp1[0], &derr[is * params.ng * params.nt + it * params.ng], gxz, params.ng, params.nz, true);
      step_forward(gp0[0], gp1[0], gp2[0], vv[0], dtz, dtx, params.nz, params.nx);

      cal_gradient(g1[0], lap[0], gp1[0], params.nz, params.nx);
      ptr = sp0;
      sp0 = sp1;
      sp1 = sp2;
      sp2 = ptr;

      ptr = gp0;
      gp0 = gp1;
      gp1 = gp2;
      gp2 = ptr;

//        float sum_g1 = std::accumulate(g1[0], g1[0] + params.nz * params.nx, 0.0f);
//        DEBUG() << format("shot %d, nt %d, sum_g1: %.6f") % (ik * size + rank) % it % sum_g1;
    }
  } /// output: derr, g1, illum
}



int main(int argc, char *argv[]) {
  /* initialize Madagascar */
  sf_init(argc, argv);

  Logger::instance().init("serial-fwi");

  DEBUG() << "HELLO";

  GlobalParams &params = GlobalParams::instance();

  float dtx = params.dt / params.dx;
  float dtz = params.dt / params.dz;
  bool csdgather = (params.csd > 0) ? true : false;

  float **vv = sf_floatalloc2(params.nz, params.nx); /* updated velocity */
  float **vtmp = sf_floatalloc2(params.nz, params.nx); /* temporary velocity computed with epsil */
  float **sp0 = sf_floatalloc2(params.nz, params.nx); /* source wavefield p0 */
  float **sp1 = sf_floatalloc2(params.nz, params.nx); /* source wavefield p1 */
  float **sp2 = sf_floatalloc2(params.nz, params.nx); /* source wavefield p2 */
  float **gp0 = sf_floatalloc2(params.nz, params.nx); /* geophone/receiver wavefield p0 */
  float **gp1 = sf_floatalloc2(params.nz, params.nx); /* geophone/receiver wavefield p1 */
  float **gp2 = sf_floatalloc2(params.nz, params.nx); /* geophone/receiver wavefield p2 */
  float **g0 = sf_floatalloc2(params.nz, params.nx); /* gradient at previous step */
  float **g1 = sf_floatalloc2(params.nz, params.nx); /* gradient at curret step */
  float **cg = sf_floatalloc2(params.nz, params.nx); /* conjugate gradient */
  float **lap = sf_floatalloc2(params.nz, params.nx); /* laplace of the source wavefield */
  float **illum = sf_floatalloc2(params.nz, params.nx); /* illumination of the source wavefield */
  float *trans = (float *)malloc(params.ng * params.nt * sizeof(float)); /* transposed one shot */

  /* initialize varibles */
  sf_floatread(vv[0], params.nz * params.nx, params.vinit);
  memset(sp0[0], 0, params.nz * params.nx * sizeof(float));
  memset(sp1[0], 0, params.nz * params.nx * sizeof(float));
  memset(sp2[0], 0, params.nz * params.nx * sizeof(float));
  memset(gp0[0], 0, params.nz * params.nx * sizeof(float));
  memset(gp1[0], 0, params.nz * params.nx * sizeof(float));
  memset(gp2[0], 0, params.nz * params.nx * sizeof(float));
  memset(g0[0], 0, params.nz * params.nx * sizeof(float));
  memset(g1[0], 0, params.nz * params.nx * sizeof(float));
  memset(cg[0], 0, params.nz * params.nx * sizeof(float));
  memset(lap[0], 0, params.nz * params.nx * sizeof(float));
  memset(vtmp[0], 0, params.nz * params.nx * sizeof(float));
  memset(illum[0], 0, params.nz * params.nx * sizeof(float));

  float *wlt = (float *)malloc(params.nt * sizeof(float)); /* ricker wavelet */
  for (int it = 0; it < params.nt; it++) {
    float tmp = SF_PI * params.fm * (it * params.dt - 1.0 / params.fm);
    tmp *= tmp;
    wlt[it] = (1.0 - 2.0 * tmp) * expf(-tmp);
  }

  if (!(params.sxbeg >= 0 && params.szbeg >= 0 && params.sxbeg + (params.ns - 1)*params.jsx < params.nx && params.szbeg + (params.ns - 1)*params.jsz < params.nz)) {
    sf_warning("sources exceeds the computing zone!\n");
    exit(1);
  }

  int *sxz = (int *)malloc(params.ns * sizeof(int)); /* source positions */
  sg_init(sxz, params.szbeg, params.sxbeg, params.jsz, params.jsx, params.ns, params.nz);

  int distx = params.sxbeg - params.gxbeg;
  int distz = params.szbeg - params.gzbeg;
  if (csdgather)  {
    if (!(params.gxbeg >= 0 && params.gzbeg >= 0 && params.gxbeg + (params.ng - 1)*params.jgx < params.nx && params.gzbeg + (params.ng - 1)*params.jgz < params.nz &&
        (params.sxbeg + (params.ns - 1)*params.jsx) + (params.ng - 1)*params.jgx - distx < params.nx  && (params.szbeg + (params.ns - 1)*params.jsz) + (params.ng - 1)*params.jgz - distz < params.nz)) {
      sf_warning("geophones exceeds the computing zone!\n");
      exit(1);
    }
  } else {
    if (!(params.gxbeg >= 0 && params.gzbeg >= 0 && params.gxbeg + (params.ng - 1)*params.jgx < params.nx && params.gzbeg + (params.ng - 1)*params.jgz < params.nz)) {
      sf_warning("geophones exceeds the computing zone!\n");
      exit(1);
    }
  }

  int *gxz = (int *)malloc(params.ng * sizeof(int)); /* geophone positions */
  sg_init(gxz, params.gzbeg, params.gxbeg, params.jgz, params.jgx, params.ng, params.nz);

  float *bndr = (float *)malloc(params.nt * (2 * params.nz + params.nx) * sizeof(float)); /* boundaries for wavefield reconstruction */
  memset(bndr, 0, params.nt * (2 * params.nz + params.nx)*sizeof(float));

  float *dobs = (float *)malloc(params.ng * params.nt * sizeof(float)); /* observed seismic data */
  memset(dobs, 0, params.ng * params.nt * sizeof(float));

  float *dcal = (float *)malloc(params.ng * sizeof(float)); /* calculated/synthetic seismic data */
  memset(dcal, 0, params.ng * sizeof(float));

  float *derr = (float *)malloc(params.ns * params.ng * params.nt * sizeof(float)); /* residual/error between synthetic and observation */
  memset(derr, 0, params.ns * params.ng * params.nt * sizeof(float));

  float *alpha1 = (float *)malloc(params.ng * sizeof(float)); /* numerator of alpha, length=ng */
  memset(alpha1, 0, params.ng * sizeof(float));

  float *alpha2 = (float *)malloc(params.ng * sizeof(float)); /* denominator of alpha, length=ng */
  memset(alpha2, 0, params.ng * sizeof(float));

  float *objval = (float *)malloc(params.niter * sizeof(float)); /* objective/misfit function */
  memset(objval, 0, params.niter * sizeof(float));

  float obj1 = 0;

  DEBUG() << "csdgather: " << csdgather;
  for (int iter = 0; iter < params.niter; iter++) {
    boost::timer::auto_cpu_timer t;
    sf_seek(params.shots, 0L, SEEK_SET);
    memcpy(g0[0], g1[0], params.nz * params.nx * sizeof(float));
    memset(g1[0], 0, params.nz * params.nx * sizeof(float));
    memset(illum[0], 0, params.nz * params.nx * sizeof(float));
//    for (int is = 0; is < params.ns; is++) {
//      sf_floatread(trans, params.ng * params.nt, params.shots);
//      matrix_transpose(trans, dobs, params.nt, params.ng);
//
////      float sum_dobs = std::accumulate(dobs, dobs + params.nt * params.ng, 0.0f);
////      DEBUG() << format("iter %d, shot %d, sum_obs %.6f") % iter % is % sum_dobs;
//
//      if (csdgather)  {
//        params.gxbeg = params.sxbeg + is * params.jsx - distx;
//        sg_init(gxz, params.gzbeg, params.gxbeg, params.jgz, params.jgx, params.ng, params.nz);
//      }
//      memset(sp0[0], 0, params.nz * params.nx * sizeof(float));
//      memset(sp1[0], 0, params.nz * params.nx * sizeof(float));
//      for (int it = 0; it < params.nt; it++) {
//        add_source(sp1, &wlt[it], &sxz[is], 1, params.nz, true);
//        step_forward(sp0, sp1, sp2, vv, dtz, dtx, params.nz, params.nx);
//
//        float **ptr = sp0;
//        sp0 = sp1;
//        sp1 = sp2;
//        sp2 = ptr;
//
//        rw_bndr(&bndr[it * (2 * params.nz + params.nx)], sp0, params.nz, params.nx, true);
//
//        record_seis(dcal, gxz, sp0, params.ng, params.nz);
//        cal_residuals(dcal, &dobs[it * params.ng], &derr[is * params.ng * params.nt + it * params.ng], params.ng);
//      }
//
//      float sum_derr = std::accumulate(&derr[is * params.ng * params.nt],
//          &derr[is * params.ng * params.nt] + params.ng * params.nt, 0.0f);
//      DEBUG() << format("shot: %d, sum_derr %.6f") % is % sum_derr;
//
//      float **ptr = sp0;
//      sp0 = sp1;
//      sp1 = ptr;
//      memset(gp0[0], 0, params.nz * params.nx * sizeof(float));
//      memset(gp1[0], 0, params.nz * params.nx * sizeof(float));
//      for (int it = params.nt - 1; it > -1; it--) {
//        rw_bndr(&bndr[it * (2 * params.nz + params.nx)], sp1, params.nz, params.nx, false);
//        step_backward(illum, lap, sp0, sp1, sp2, vv, dtz, dtx, params.nz, params.nx);
//        add_source(sp1, &wlt[it], &sxz[is], 1, params.nz, false);
//
//        add_source(gp1, &derr[is * params.ng * params.nt + it * params.ng], gxz, params.ng, params.nz, true);
//        step_forward(gp0, gp1, gp2, vv, dtz, dtx, params.nz, params.nx);
//
//        cal_gradient(g1, lap, gp1, params.nz, params.nx);
//        ptr = sp0;
//        sp0 = sp1;
//        sp1 = sp2;
//        sp2 = ptr;
//        ptr = gp0;
//        gp0 = gp1;
//        gp1 = gp2;
//        gp2 = ptr;
//
////        float sum_g1 = std::accumulate(g1[0], g1[0] + params.nz * params.nx, 0.0f);
////        DEBUG() << format("shot %d, nt %d, sum_g1: %.6f") % (is) % it % sum_g1;
//      }
//    }
    {
      float sum_vv = std::accumulate(vv[0], vv[0] + params.nx * params.nz, 0.0f);
      DEBUG() << format("iter %d, sum_vv %.6f") % iter % sum_vv;
    }
    cal_derr_illum_grad(params, vv, derr, illum, g1, wlt, sxz, gxz);
    float obj = cal_objective(derr, params.ng * params.nt * params.ns);

    float sum_g1 = std::accumulate(g1[0], g1[0] + params.nz * params.nx, 0.0f);
    DEBUG() << format("iter: %d, sum_g1: %.6f") % iter % sum_g1;

    float sum_illum = std::accumulate(illum[0], illum[0] + params.nz * params.nx, 0.0f);
    DEBUG() << format("iter: %d, sum_illum: %.6f") % iter % sum_illum;

    scale_gradient(g1[0], vv[0], illum[0], params.nz, params.nx, params.precon);
    sf_floatwrite(illum[0], params.nz * params.nx, params.illums);
    bell_smoothz(g1[0], illum[0], params.rbell, params.nz, params.nx);
    bell_smoothx(illum[0], g1[0], params.rbell, params.nz, params.nx);
    sf_floatwrite(g1[0], params.nz * params.nx, params.grads);


    float beta;
    if (iter > 0) {
      beta = cal_beta(g0[0], g1[0], cg[0], params.nz, params.nx);
    } else {
      beta = 0.0;
    }
    cal_conjgrad(g1[0], cg[0], beta, params.nz, params.nx);
    float epsil = cal_epsilon(vv[0], cg[0], params.nz, params.nx);

    sf_seek(params.shots, 0L, SEEK_SET);
    memset(alpha1, 0, params.ng * sizeof(float));
    memset(alpha2, 0, params.ng * sizeof(float));
    cal_vtmp(vtmp[0], vv[0], cg[0], epsil, params.nz, params.nx);

    float sum_vtmp = std::accumulate(vtmp[0], vtmp[0] + params.nz * params.nx, 0.0f);
    DEBUG() << format("iter %d, sum_vtmp %.0f") % iter % sum_vtmp;

    for (int is = 0; is < params.ns; is++) {
      sf_floatread(trans, params.ng * params.nt, params.shots);
      matrix_transpose(trans, dobs, params.nt, params.ng);

//      float sum_dobs = std::accumulate(dobs, dobs + params.nt * params.ng, 0.0f);
//      DEBUG() << format("iter %d, shot %d, sum_obs %.6f") % iter % is % sum_dobs;

      if (csdgather)  {
        params.gxbeg = params.sxbeg + is * params.jsx - distx;
        sg_init(gxz, params.gzbeg, params.gxbeg, params.jgz, params.jgx, params.ng, params.nz);
      }
      memset(sp0[0], 0, params.nz * params.nx * sizeof(float));
      memset(sp1[0], 0, params.nz * params.nx * sizeof(float));
      for (int it = 0; it < params.nt; it++) {
        add_source(sp1[0], &wlt[it], &sxz[is], 1, params.nz, true);
        step_forward(sp0[0], sp1[0], sp2[0], vtmp[0], dtz, dtx, params.nz, params.nx);
        float **ptr = sp0;
        sp0 = sp1;
        sp1 = sp2;
        sp2 = ptr;

        record_seis(dcal, gxz, sp0[0], params.ng, params.nz);
        sum_alpha12(alpha1, alpha2, dcal, &dobs[it * params.ng], &derr[is * params.ng * params.nt + it * params.ng], params.ng);
      }
    }

    float sum_alpha1 = std::accumulate(alpha1, alpha1 + params.ng, 0.0f);
    DEBUG() << format("iter %d, sum_alpha1 %.6f") % iter % sum_alpha1;
    float sum_alpha2 = std::accumulate(alpha2, alpha2 + params.ng, 0.0f);
    DEBUG() << format("iter %d, sum_alpha2 %.6f") % iter % sum_alpha2;

    float alpha = cal_alpha(alpha1, alpha2, epsil, params.ng);
    update_vel(vv[0], cg[0], alpha, params.nz, params.nx);

    float sum_vv = std::accumulate(vv[0], vv[0] + params.nx * params.nz, 0.0f);
    DEBUG() << format("iter %d, sum_vv %.6f") % iter % sum_vv;

    sf_floatwrite(vv[0], params.nz * params.nx, params.vupdates);

    if (iter == 0) {
      obj1 = obj;
      objval[iter] = 1.0;
    } else {
      objval[iter] = obj / obj1;
    }

    if (params.verb) { // output important information at each FWI iteration
      INFO() << format("iteration %d obj=%f  beta=%f  epsil=%f  alpha=%f") % (iter + 1) % obj % beta % epsil % alpha;
    }

  }
  sf_floatwrite(objval, params.niter, params.objs);

  sf_close();
  free(*vv);
  free(vv);
  free(*vtmp);
  free(vtmp);
  free(*sp0);
  free(sp0);
  free(*sp1);
  free(sp1);
  free(*sp2);
  free(sp2);
  free(*gp0);
  free(gp0);
  free(*gp1);
  free(gp1);
  free(*gp2);
  free(gp2);
  free(*g0);
  free(g0);
  free(*g1);
  free(g1);
  free(*cg);
  free(cg);
  free(*lap);
  free(lap);
  free(*illum);
  free(illum);
  free(objval);
  free(wlt);
  free(sxz);
  free(gxz);
  free(bndr);
  free(trans);
  free(dobs);
  free(dcal);
  free(derr);
  free(alpha1);
  free(alpha2);

  exit(0);
}
