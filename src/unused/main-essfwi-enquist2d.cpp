
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
#include "essfwi-params.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "cycle-swap.h"
#include "sum.h"
#include "sf-velocity-reader.h"
#include "shotdata-reader.h"
#include "velocity.h"
#include "shot-position.h"
#include "random-code.h"
#include "encoder.h"
#include "enquistabc2d.h"

static float cal_obj_derr_illum_grad(const EssFwiParams &params,
    std::vector<float> &derr,  /* output */
    std::vector<float> &illum, /* output */
    std::vector<float> &g1,    /* output */
    const std::vector<float> &encsrc,
    const std::vector<float> &encobs,
    const EnquistAbc2d &fmMethod,
    const ShotPosition &allSrcPos,
    const ShotPosition &allGeoPos)
{
  int nt = params.nt;
  int nz = params.nz;
  int nx = params.nx;
  int ng = params.ng;
  int ns = params.ns;

  std::vector<float> bndr = fmMethod.initBndryVector(nt);
  std::vector<float> dcal(ng, 0); /* calculated/synthetic seismic data */

  std::vector<float> sp0(nz * nx); /* source wavefield p0 */
  std::vector<float> sp1(nz * nx); /* source wavefield p1 */
  std::vector<float> sp2(nz * nx); /* source wavefield p2 */
  std::vector<float> gp0(nz * nx); /* geophone/receiver wavefield p0 */
  std::vector<float> gp1(nz * nx); /* geophone/receiver wavefield p1 */
  std::vector<float> gp2(nz * nx); /* geophone/receiver wavefield p2 */
  std::vector<float> lap(nz * nx); /* laplace of the source wavefield */


  for (int it = 0; it < nt; it++) {
    fmMethod.addSource(&sp1[0], &encsrc[it * ns], allSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0], &sp2[0]);

    // cycle swap
    cycleSwap(sp0, sp1, sp2);

    fmMethod.writeBndry(&bndr[0], &sp0[0], it);
    fmMethod.recordSeis(&dcal[0], &sp0[0], allGeoPos);

    cal_residuals(&dcal[0], &encobs[it * ng], &derr[it * ng], ng);
  }

  std::swap(sp0, sp1);


  for (int it = nt - 1; it > -1; it--) {
    /// backward propagate source wavefield
    fmMethod.readBndry(&bndr[0], &sp1[0], it);
    fmMethod.stepBackward(&illum[0], &lap[0], &sp0[0], &sp1[0], &sp2[0]);
    fmMethod.subSource(&sp1[0], &encsrc[it * ns], allSrcPos);

    /// forward propagate geophone wavefield
    fmMethod.addSource(&gp1[0], &derr[it * ng], allGeoPos);
    fmMethod.stepForward(&gp0[0], &gp1[0], &gp2[0]);

    /// calculate gradient
    cal_gradient(&g1[0], &lap[0], &gp1[0], nz, nx);

    cycleSwap(sp0, sp1, sp2);
    cycleSwap(gp0, gp1, gp2);
  }


  float obj = cal_objective(&derr[0], ng * nt);

  return obj;
}

float calVelUpdateStepLen(const EssFwiParams &params,
    const std::vector<float> &encsrc,
    const std::vector<float> &encobs,
    const std::vector<float> &derr,
    float epsil,
    const EnquistAbc2d &fmMethod,
    const ShotPosition &allSrcPos,
    const ShotPosition &allGeoPos
    )
{
  int nt = params.nt;
  int nz = params.nz;
  int nx = params.nx;
  int ng = params.ng;
  int ns = params.ns;

  std::vector<float> dcal(ng, 0); /* calculated/synthetic seismic data */
  std::vector<float> sp0(nz * nx); /* source wavefield p0 */
  std::vector<float> sp1(nz * nx); /* source wavefield p1 */
  std::vector<float> sp2(nz * nx); /* source wavefield p2 */

  std::vector<float> alpha1(ng, 0); /* numerator of alpha, length=ng */
  std::vector<float> alpha2(ng, 0); /* denominator of alpha, length=ng */

  for (int it = 0; it < nt; it++) {
    fmMethod.addSource(&sp1[0], &encsrc[it * ns], allSrcPos);
    fmMethod.stepForward(&sp0[0], &sp1[0], &sp2[0]);

    cycleSwap(sp0, sp1, sp2);

    fmMethod.recordSeis(&dcal[0], &sp0[0], allGeoPos);

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

  int ns = params.ns;
  int nt = params.nt;
  int ng = params.ng;

  std::vector<float> dobs(params.ns * params.nt * params.ng); /* observed data */
  std::vector<float> cg(params.nz * params.nx, 0);    /* conjugate gradient */
  std::vector<float> g0(params.nz * params.nx, 0);    /* gradient at previous step */
  std::vector<float> wlt(params.nt); /* ricker wavelet */
  std::vector<float> objval(params.niter, 0); /* objective/misfit function */

  /* initialize wavelet */
  rickerWavelet(&wlt[0], params.nt, params.fm, params.dt, params.amp);

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, params.ns, params.nz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, params.ng, params.nz);

  // read velocity
  Velocity v0 = SfVelocityReader::read(params.vinit, params.nx, params.nz);

  // read observed data
  ShotDataReader::serialRead(params.shots, &dobs[0], params.ns, params.nt, params.ng);

  EnquistAbc2d fmMethod(params.dt, params.dx, params.dz);
  Velocity vel = fmMethod.expandDomain(v0);

  float obj0 = 0;
  for (int iter = 0; iter < params.niter; iter++) {
    boost::timer::cpu_timer timer;
    std::vector<float> g1(params.nz * params.nx, 0);    /* gradient at curret step */
    std::vector<float> derr(params.ng * params.nt, 0); /* residual/error between synthetic and observation */
    std::vector<float> illum(params.nz * params.nx, 0); /* illumination of the source wavefield */
    Velocity vtmp = vel;  /* temporary velocity computed with epsil */

    fmMethod.bindVelocity(vel);

    // create random codes
    const std::vector<int> encodes = RandomCode::genPlus1Minus1(ns);
    std::copy(encodes.begin(), encodes.end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

    Encoder encoder(encodes);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    /**
     * calculate local objective function & derr & illum & g1(gradient)
     */
    float obj = cal_obj_derr_illum_grad(params, derr, illum, g1, encsrc, encobs, fmMethod, allSrcPos, allGeoPos);

//    DEBUG() << format("sum_derr %f, sum_illum %f, sum_g1 %f") % sum(derr) % sum(illum) % sum(g1);
    objval[iter] = iter == 0 ? obj0 = obj, 1.0 : obj / obj0;

    float epsil = 0;
    float beta = 0;
    sf_floatwrite(&illum[0], params.nz * params.nx, params.illums);

    scale_gradient(&g1[0], &vel.dat[0], &illum[0], params.nz, params.nx, params.precon);
    bell_smoothz(&g1[0], &illum[0], params.rbell, params.nz, params.nx);
    bell_smoothx(&illum[0], &g1[0], params.rbell, params.nz, params.nx);
    sf_floatwrite(&g1[0], params.nz * params.nx, params.grads);

//    DEBUG() << format("before beta: sum_g0: %f, sum_g1: %f, sum_cg: %f") % sum(g0) % sum(g1) % sum(cg);
    beta = iter == 0 ? 0.0 : cal_beta(&g0[0], &g1[0], &cg[0], params.nz, params.nx);

    cal_conjgrad(&g1[0], &cg[0], beta, params.nz, params.nx);
    epsil = cal_epsilon(&vel.dat[0], &cg[0], params.nz, params.nx);
    cal_vtmp(&vtmp.dat[0], &vel.dat[0], &cg[0], epsil, params.nz, params.nx);

    std::swap(g1, g0); // let g0 be the previous gradient

    fmMethod.bindVelocity(vtmp);
    float alpha = calVelUpdateStepLen(params, encsrc, encobs, derr, epsil, fmMethod, allSrcPos, allGeoPos);

    update_vel(&vel.dat[0], &cg[0], alpha, params.nz, params.nx);

    sf_floatwrite(&vel.dat[0], params.nz * params.nx, params.vupdates);

    // output important information at each FWI iteration
    INFO() << format("iteration %d obj=%e  beta=%f  epsil=%f  alpha=%f") % (iter + 1) % obj % beta % epsil % alpha;
//    INFO() << timer.format(2);

  } /// end of iteration

  sf_floatwrite(&objval[0], params.niter, params.objs);

  sf_close();

  return 0;
}
