
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
#include <numeric>
#include <cstdlib>
#include <vector>

#include <boost/timer/timer.hpp>
#include "logger.h"
#include "fwi-params.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "cycle-swap.h"
#include "sum.h"
#include "sf-velocity-reader.h"
#include "shotdata-reader.h"
#include "velocity.h"
#include "shot-position.h"
#include "spongeabc4d.h"

float cal_obj_derr_illum_grad(
    float *derr,  /* output */
    float *illum, /* output */
    float *g1,    /* output */
    const float *wlt,
    const float *dobs,
    int nt,
    const SpongeAbc4d &fmMethod,
    const ShotPosition &allSrcPos,
    const ShotPosition &allGeoPos)
{
  int nz = fmMethod.getVelocity().nz;
  int nx = fmMethod.getVelocity().nx;
  int ng = allGeoPos.ns;
  int ns = allSrcPos.ns;

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> dcal(ng, 0); /* calculated/synthetic seismic data */
  std::vector<float> sp0(nz * nx, 0); /* source wavefield p0 */
  std::vector<float> sp1(nz * nx, 0); /* source wavefield p1 */
  std::vector<float> gp0(nz * nx, 0); /* geophone/receiver wavefield p0 */
  std::vector<float> gp1(nz * nx, 0); /* geophone/receiver wavefield p1 */
  std::vector<float> lap(nz * nx, 0); /* laplace of the source wavefield */

  for (int is = 0; is < ns; is++) {
    std::fill(sp0.begin(), sp0.end(), 0);
    std::fill(sp1.begin(), sp1.end(), 0);

    ShotPosition curSrcPos = allSrcPos.clip(is);

    for (int it = 0; it < nt; it++) {
      fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
      fmMethod.stepForward(&sp0[0], &sp1[0]);
      fmMethod.applySponge(&sp0[0]);
      fmMethod.applySponge(&sp1[0]);

      swap(sp0, sp1);

      fmMethod.writeBndry(&bndr[0], &sp0[0], it);
      fmMethod.recordSeis(&dcal[0], &sp0[0], allGeoPos);

      cal_residuals(&dcal[0], &dobs[is * nt * ng + it * ng], &derr[is * ng * nt + it * ng], ng);
    }

    std::fill(gp0.begin(), gp0.end(), 0);
    std::fill(gp1.begin(), gp1.end(), 0);

    for (int it = nt - 1; it > -1; it--) {
      /// backward propagate source wavefield
      fmMethod.readBndry(&bndr[0], &sp0[0], it);
      std::swap(sp0, sp1);
      fmMethod.stepBackward(illum, &lap[0], &sp0[0], &sp1[0]);
      fmMethod.subSource(&sp1[0], &wlt[it], curSrcPos);

      /// forward propagate geophone wavefield
      fmMethod.addSource(&gp1[0], &derr[is * ng * nt + it * ng], allGeoPos);
      fmMethod.stepForward(&gp0[0], &gp1[0]);
      fmMethod.applySponge(&gp0[0]);
      fmMethod.applySponge(&gp1[0]);

      /// calculate gradient
      cal_gradient(&g1[0], &lap[0], &gp1[0], nz, nx);

      std::swap(gp0, gp1);
    }

  } /// output: derr, g1, illum

  float obj = cal_objective(&derr[0], ng * nt * ns);

  return obj;
}

float calVelUpdateStepLen(
    const float *wlt,
    const float *dobs,
    const float *derr,
    float epsil,
    int nt,
    const SpongeAbc4d &fmMethod,
    const ShotPosition &allSrcPos,
    const ShotPosition &allGeoPos
    )
{
  int nz = fmMethod.getVelocity().nz;
  int nx = fmMethod.getVelocity().nx;
  int ng = allGeoPos.ns;
  int ns = allSrcPos.ns;

  std::vector<float> dcal(ng, 0); /* calculated/synthetic seismic data */
  std::vector<float> sp0(nz * nx, 0); /* source wavefield p0 */
  std::vector<float> sp1(nz * nx, 0); /* source wavefield p1 */

  std::vector<float> alpha1(ng, 0); /* numerator of alpha, length=ng */
  std::vector<float> alpha2(ng, 0); /* denominator of alpha, length=ng */

  for (int is = 0; is < ns; is++) {
    std::fill(sp0.begin(), sp0.end(), 0);
    std::fill(sp1.begin(), sp1.end(), 0);
    ShotPosition curSrcPos = allSrcPos.clip(is);

    for (int it = 0; it < nt; it++) {
      fmMethod.addSource(&sp1[0], &wlt[it], curSrcPos);
      fmMethod.stepForward(&sp0[0], &sp1[0]);
      fmMethod.applySponge(&sp0[0]);
      fmMethod.applySponge(&sp1[0]);

      swap(sp0, sp1);

      fmMethod.recordSeis(&dcal[0], &sp0[0], allGeoPos);

      sum_alpha12(&alpha1[0], &alpha2[0], &dcal[0], &dobs[is * nt * ng + it * ng], &derr[is * ng * nt + it * ng], ng);
    }
  }

  float alpha = cal_alpha(&alpha1[0], &alpha2[0], epsil, ng);

  return alpha;
}

int main(int argc, char *argv[]) {

  /* initialize Madagascar */
  sf_init(argc, argv);

  Logger::instance().init("serial-fwi");

  FwiParams &params = FwiParams::instance();

  int ns = params.ns;   // # of shots
  int ng = params.ng;   // # of geophones
  int nt = params.nt;   // # of time steps
  int niter = params.niter; // # of iteration
  int mnz = params.nz;  // model nz;
  int mnx = params.nx;  // model nx
  int nb = params.nb;   // # layer of boundary
  float dt = params.dt;
  float dx = params.dx;
  float dz = params.dz;

  DEBUG() << "nb: " << params.nb;

  /* initialize wavelet */
  std::vector<float> wlt(nt); /* ricker wavelet */
  rickerWavelet(&wlt[0], nt, params.fm, dt, params.amp);

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, mnz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, mnz);

  // read velocity
  Velocity v0 = SfVelocityReader::read(params.vinit, mnx, mnz);

  // read observed data
  std::vector<float> dobs(ns * nt * ng); /* observed data */
  ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);

  SpongeAbc4d fmMethod(dt, dx, dz, nb);
  Velocity exvel = fmMethod.expandDomain(v0);

  int nxpad = exvel.nx;
  int nzpad = exvel.nz;
  std::vector<float> cg(nxpad * nzpad, 0);    /* conjugate gradient */
  std::vector<float> g0(nxpad * nzpad, 0);    /* gradient at previous step */
  std::vector<float> objval(niter, 0);        /* objective/misfit function */

  float obj0 = 0;
  for (int iter = 0; iter < niter; iter++) {
    boost::timer::cpu_timer timer;
    std::vector<float> g1(nxpad * nzpad, 0);    /* gradient at curret step */
    std::vector<float> derr(ns * ng * nt, 0);   /* residual/error between synthetic and observation */
    std::vector<float> illum(nxpad * nxpad, 0); /* illumination of the source wavefield */
    std::vector<float> tmpForWrite(mnx * mnz);  /* model size for write */
    Velocity vtmp = exvel;                      /* temporary velocity computed with epsil */
    fmMethod.bindVelocity(exvel);

    /**
     * calculate local objective function & derr & illum & g1(gradient)
     */
    float obj = cal_obj_derr_illum_grad(&derr[0], &illum[0], &g1[0], &wlt[0], &dobs[0], params.nt, fmMethod, allSrcPos, allGeoPos);

    DEBUG() << format("sum_derr %f, sum_illum %f, sum_g1 %f") % sum(derr) % sum(illum) % sum(g1);

    objval[iter] = iter == 0 ? obj0 = obj, 1.0 : obj / obj0;

    float epsil = 0;
    float beta = 0;
    fmMethod.shrinkDomain(&tmpForWrite[0], &illum[0], mnx, mnz);
    sf_floatwrite(&tmpForWrite[0], tmpForWrite.size(), params.illums);

    scale_gradient(&g1[0], &exvel.dat[0], &illum[0], nzpad, nxpad, params.precon);
    bell_smoothz(&g1[0], &illum[0], params.rbell, nzpad, nxpad);
    bell_smoothx(&illum[0], &g1[0], params.rbell, nzpad, nxpad);

    fmMethod.shrinkDomain(&tmpForWrite[0], &g1[0], mnx, mnz);
    sf_floatwrite(&tmpForWrite[0], tmpForWrite.size(), params.grads);

    DEBUG() << format("before beta: sum_g0: %f, sum_g1: %f, sum_cg: %f") % sum(g0) % sum(g1) % sum(cg);
    beta = iter == 0 ? 0.0 : cal_beta(&g0[0], &g1[0], &cg[0], nzpad, nxpad);

    cal_conjgrad(&g1[0], &cg[0], beta, nzpad, nxpad);
    epsil = cal_epsilon(&exvel.dat[0], &cg[0], nzpad, nxpad);
    cal_vtmp(&vtmp.dat[0], &exvel.dat[0], &cg[0], epsil, nzpad, nxpad);

    std::swap(g1, g0); // let g0 be the previous gradient

    fmMethod.bindVelocity(vtmp);
    float alpha = calVelUpdateStepLen(&wlt[0], &dobs[0], &derr[0], epsil, nt, fmMethod, allSrcPos, allGeoPos);

    update_vel(&exvel.dat[0], &cg[0], alpha, nzpad, nxpad);

    fmMethod.shrinkDomain(&tmpForWrite[0], &exvel.dat[0], mnx, mnz);
    sf_floatwrite(&tmpForWrite[0], tmpForWrite.size(), params.vupdates);

    // output important information at each FWI iteration
    INFO() << format("iteration %d obj=%f  beta=%f  epsil=%f  alpha=%f") % (iter + 1) % obj % beta % epsil % alpha;
    INFO() << timer.format(2);

  } /// end of iteration

  sf_floatwrite(&objval[0], params.niter, params.objs);

  sf_close();

  return 0;
}
