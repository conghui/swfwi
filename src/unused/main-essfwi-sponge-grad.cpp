
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
#include "zjh4t10dspongevelnotrans.h"
#include "sfutil.h"
#include "aux.h"
#include "preserved-alpha.h"


static float cal_obj_derr_illum_grad(
    std::vector<float> &derr,  /* output */
    std::vector<float> &illum, /* output */
    std::vector<float> &g1,    /* output */
    const std::vector<float> &encsrc,
    const std::vector<float> &encobs,
    int nt,
    const Zjh4t10dSpongeVelNoTrans &fmMethod,
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

  for (int it = 0; it < nt; it++) {
    fmMethod.addSource(&sp1[0], &encsrc[it * ns], allSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);

    swap(sp0, sp1);

    fmMethod.recordSeis(&dcal[0], &sp0[0]);
    cal_residuals(&dcal[0], &encobs[it * ng], &derr[it * ng], ng);
    fmMethod.writeBndry(&bndr[0], &sp0[0], it);
  }

  for (int it = nt - 1; it > -1; it--) {
    /// backward propagate source wavefield
    fmMethod.readBndry(&bndr[0], &sp0[0], it);
    std::swap(sp0, sp1);
    fmMethod.stepBackwardLapIllum(&lap[0], &illum[0], &sp0[0], &sp1[0]);
    fmMethod.subEncodedSource(&sp0[0], &encsrc[it * ns]);

    /// forward propagate geophone wavefield
    fmMethod.addSource(&gp1[0], &derr[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0]);

    /// calculate gradient
    cal_gradient(&g1[0], &lap[0], &gp1[0], nz, nx);

    std::swap(gp0, gp1);
  }


  float obj = cal_objective(&derr[0], ng * nt);

  return obj;
}

float calVelUpdateStepLen(
    const std::vector<float> &encsrc,
    const std::vector<float> &encobs,
    const std::vector<float> &derr,
    float epsil,
    int nt,
    const Zjh4t10dSpongeVelNoTrans &fmMethod,
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


  for (int it = 0; it < nt; it++) {
    fmMethod.addSource(&sp1[0], &encsrc[it * ns], allSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0]);

    swap(sp0, sp1);

    fmMethod.recordSeis(&dcal[0], &sp0[0]);
    sum_alpha12(&alpha1[0], &alpha2[0], &dcal[0], &encobs[it * ng], &derr[it * ng], ng);
  }

  float alpha = cal_alpha(&alpha1[0], &alpha2[0], epsil, ng);

  return alpha;
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

  Zjh4t10dSpongeVelNoTrans fmMethod(allSrcPos, allGeoPos, dt, dx, fm, nb, nt);

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);
  Velocity exvel = fmMethod.expandDomain(v0);

  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);


  int nzpad = exvel.nz;
  int nxpad = exvel.nx;
  std::vector<float> dobs(ns * nt * ng);     /* all observed data */
  std::vector<float> g0(exvel.nx * exvel.nz, 0); /* gradient at previous step */

  ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);
  //sfFloatWrite1d("orgdata.rsf", &dobs[0], ns * nt * ng);

  std::vector<float> updateDirection(exvel.nx * exvel.nz, 0);
  std::vector<float> cg(nxpad * nzpad, 0);    /* conjugate gradient */
  std::vector<float> objval(nt, 0);
  float obj0 = 0;

  for (int iter = 0; iter < params.niter; iter++) {
    boost::timer::cpu_timer timer;

    std::vector<float> derr(ng * nt, 0);   /* residual/error between synthetic and observation */
    std::vector<float> illum(exvel.nz * exvel.nx, 0); /* illumination of the source wavefield */
    std::vector<float> g1(exvel.nx * exvel.nz, 0); /* gradient at current step */
    Velocity vtmp = exvel;                      /* temporary velocity computed with epsil */

    // create random codes
    const std::vector<int> encodes = RandomCode::genPlus1Minus1(params.ns);
//    std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

    Encoder encoder(encodes);
    std::vector<float> encobs = encoder.encodeObsData(dobs, params.nt, params.ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    float obj = cal_obj_derr_illum_grad(derr, illum, g1,encsrc, encobs, nt, fmMethod, allSrcPos, allGeoPos);

    objval[iter] = iter == 0 ? obj0 = obj, 1.0 : obj / obj0;

    float epsil = 0;
    float beta = 0;
    fmMethod.sfWrite(illum, params.illums);

    scale_gradient(&g1[0], &exvel.dat[0], &illum[0], nzpad, nxpad, params.precon);
    bell_smoothz(&g1[0], &illum[0], params.rbell, nzpad, nxpad);
    bell_smoothx(&illum[0], &g1[0], params.rbell, nzpad, nxpad);
    fmMethod.sfWrite(g1, params.grads);

    beta = iter == 0 ? 0.0 : cal_beta(&g0[0], &g1[0], &cg[0], nzpad, nxpad);

    cal_conjgrad(&g1[0], &cg[0], beta, nzpad, nxpad);
    epsil = cal_epsilon(&exvel.dat[0], &cg[0], nzpad, nxpad);
    cal_vtmp(&vtmp.dat[0], &exvel.dat[0], &cg[0], epsil, nzpad, nxpad);

    std::swap(g1, g0); // let g0 be the previous gradient

    Zjh4t10dSpongeVelNoTrans tmpfm = fmMethod;
    tmpfm.bindVelocity(vtmp);
    float alpha = calVelUpdateStepLen(encsrc, encobs, derr, epsil, nt, tmpfm, allSrcPos, allGeoPos);

    update_vel(&exvel.dat[0], &cg[0], alpha, nzpad, nxpad);

    // output important information at each FWI iteration
    INFO() << format("iteration %d obj=%e  beta=%e  epsil=%e  alpha=%e") % (iter + 1) % obj % beta % epsil % alpha;
//    INFO() << timer.format(2);

    fmMethod.sfWrite(exvel.dat, params.vupdates);
  } /// end of iteration

  sf_floatwrite(&objval[0], objval.size(), params.objs);

  sf_close();

  return 0;
}
