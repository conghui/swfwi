
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
#include <cmath>

#include <omp.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include <boost/timer/timer.hpp>
#include "logger.h"
#include "global-params.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "cycle-swap.h"

void MpiInplaceReduce(void *buf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, buf, count, datatype, op, root, comm);
  } else {
    MPI_Reduce(buf, NULL, count, datatype, op, root, comm);
  }
}


float cal_obj_derr_illum_grad(const GlobalParams &params,
    float *derr,  /* output */
    float *illum, /* output */
    float *g1,    /* output */
    const float *vv,
    const float *wlt,
    const int *sxz,
    const int *gxz)
{
  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<float> trans(params.ng * params.nt);
  std::vector<float> dobs(params.ng * params.nt);
  std::vector<float> bndr(params.nt * (2 * params.nz + params.nx), 0); /* boundaries for wavefield reconstruction */
  std::vector<float> dcal(params.ng, 0); /* calculated/synthetic seismic data */

  std::vector<float> sp0(params.nz * params.nx); /* source wavefield p0 */
  std::vector<float> sp1(params.nz * params.nx); /* source wavefield p1 */
  std::vector<float> sp2(params.nz * params.nx); /* source wavefield p2 */
  std::vector<float> gp0(params.nz * params.nx); /* geophone/receiver wavefield p0 */
  std::vector<float> gp1(params.nz * params.nx); /* geophone/receiver wavefield p1 */
  std::vector<float> gp2(params.nz * params.nx); /* geophone/receiver wavefield p2 */
  std::vector<float> lap(params.nz * params.nx); /* laplace of the source wavefield */

  float dtx = params.dt / params.dx;
  float dtz = params.dt / params.dz;

  for (int is = rank, ik = 0; is < params.ns; is += size, ik++) {
    sf_floatread(&trans[0], params.ng * params.nt, params.shots); /* Read local portion of input data */
    matrix_transpose(&trans[0], &dobs[0], params.nt, params.ng);

    std::fill(sp0.begin(), sp0.end(), 0);
    std::fill(sp1.begin(), sp1.end(), 0);
    for (int it = 0; it < params.nt; it++) {
      add_source(&sp1[0], &wlt[it], &sxz[is], 1, params.nz, true);
      step_forward(&sp0[0], &sp1[0], &sp2[0], vv, dtz, dtx, params.nz, params.nx);

      // cycle swap
      cycleSwap(sp0, sp1, sp2);

      rw_bndr(&bndr[it * (2 * params.nz + params.nx)], &sp0[0], params.nz, params.nx, true);

      record_seis(&dcal[0], gxz, &sp0[0], params.ng, params.nz);
      cal_residuals(&dcal[0], &dobs[it * params.ng], &derr[ik * params.ng * params.nt + it * params.ng], params.ng);
    }

    std::swap(sp0, sp1);

    std::fill(gp0.begin(), gp0.end(), 0);
    std::fill(gp1.begin(), gp1.end(), 0);

    for (int it = params.nt - 1; it > -1; it--) {
      rw_bndr(&bndr[it * (2 * params.nz + params.nx)], &sp1[0], params.nz, params.nx, false);
      step_backward(illum, &lap[0], &sp0[0], &sp1[0], &sp2[0], vv, dtz, dtx, params.nz, params.nx);
      add_source(&sp1[0], &wlt[it], &sxz[is], 1, params.nz, false);

      add_source(&gp1[0], &derr[ik * params.ng * params.nt + it * params.ng], gxz, params.ng, params.nz, true);
      step_forward(&gp0[0], &gp1[0], &gp2[0], vv, dtz, dtx, params.nz, params.nx);

      cal_gradient(&g1[0], &lap[0], &gp1[0], params.nz, params.nx);

      cycleSwap(sp0, sp1, sp2);
      cycleSwap(gp0, gp1, gp2);
    }

    sf_seek(params.shots, params.nt * params.ng * (size - 1)*sizeof(float), SEEK_CUR); /* Move on to the next portion */

  } /// output: derr, g1, illum

  float obj = cal_objective(&derr[0], params.ng * params.nt * params.nk);

  return obj;
}

float calVelUpdateStepLen(const GlobalParams &params,
    const float *vtmp,
    const float *wlt,
    const int *sxz,
    const int *gxz,
    const float *derr,
    float epsil
    )
{
  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<float> trans(params.ng * params.nt);
  std::vector<float> dobs(params.ng * params.nt);
  std::vector<float> dcal(params.ng, 0); /* calculated/synthetic seismic data */
  std::vector<float> sp0(params.nz * params.nx); /* source wavefield p0 */
  std::vector<float> sp1(params.nz * params.nx); /* source wavefield p1 */
  std::vector<float> sp2(params.nz * params.nx); /* source wavefield p2 */

  std::vector<float> alpha1(params.ng, 0); /* numerator of alpha, length=ng */
  std::vector<float> alpha2(params.ng, 0); /* denominator of alpha, length=ng */

  float dtx = params.dt / params.dx;
  float dtz = params.dt / params.dz;
  for (int is = rank, ik = 0; is < params.ns; is+= size, ik++) {
    sf_floatread(&trans[0], params.ng * params.nt, params.shots);
    matrix_transpose(&trans[0], &dobs[0], params.nt, params.ng);

    std::fill(sp0.begin(), sp0.end(), 0);
    std::fill(sp1.begin(), sp1.end(), 0);
    for (int it = 0; it < params.nt; it++) {
      add_source(&sp1[0], &wlt[it], &sxz[is], 1, params.nz, true);
      step_forward(&sp0[0], &sp1[0], &sp2[0], &vtmp[0], dtz, dtx, params.nz, params.nx);

      std::swap(sp0, sp1);
      std::swap(sp1, sp2);

      record_seis(&dcal[0], gxz, &sp0[0], params.ng, params.nz);
      sum_alpha12(&alpha1[0], &alpha2[0], &dcal[0], &dobs[it * params.ng], &derr[ik * params.ng * params.nt + it * params.ng], params.ng);
    }
    sf_seek(params.shots, params.nt * params.ng * (size - 1)*sizeof(float), SEEK_CUR); /* Move on to the next portion */
  }


  /**
   * MPI reduce
   */
  MpiInplaceReduce(&alpha1[0], params.ng, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MpiInplaceReduce(&alpha2[0], params.ng, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  float alpha = cal_alpha(&alpha1[0], &alpha2[0], epsil, params.ng);

  return alpha;
}



int main(int argc, char *argv[]) {
  MPI_Init (&argc, &argv);

  /* initialize Madagascar */
  sf_init(argc, argv);

  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);/* how many nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);/* who am I? */

  Logger::instance().init("mpi-fwi", boost::log::trivial::debug, Logger::NO_TIMESTAMP, rank);

  GlobalParams &params = GlobalParams::instance();

  // how many groups of MPI chunk
  INFO() << format("each process should process %d shots") % params.nk;

  std::vector<float> vv(params.nz * params.nx, 0);    /* updated velocity */
  std::vector<float> vtmp(params.nz * params.nx, 0);  /* temporary velocity computed with epsil */
  std::vector<float> g0(params.nz * params.nx, 0);    /* gradient at previous step */
  std::vector<float> g1(params.nz * params.nx, 0);    /* gradient at curret step */
  std::vector<float> cg(params.nz * params.nx, 0);    /* conjugate gradient */
  std::vector<float> illum(params.nz * params.nx, 0); /* illumination of the source wavefield */
  std::vector<float> wlt(params.nt); /* ricker wavelet */
  std::vector<int> sxz(params.ns); /* source positions */
  std::vector<int> gxz(params.ng); /* geophone positions */
  std::vector<float> objval(params.niter, 0); /* objective/misfit function */
  std::vector<float> derr(params.nk * params.ng * params.nt, 0); /* residual/error between synthetic and observation */

  /* initialize varibles */
  sf_floatread(&vv[0], params.nz * params.nx, params.vinit);
  rickerWavelet(&wlt[0], params.nt, params.fm, params.dt);
  sg_init(&sxz[0], params.szbeg, params.sxbeg, params.jsz, params.jsx, params.ns, params.nz);
  sg_init(&gxz[0], params.gzbeg, params.gxbeg, params.jgz, params.jgx, params.ng, params.nz);

  float obj0 = 0;
  for (int iter = 0; iter < params.niter; iter++) {
    boost::timer::cpu_timer timer;
    sf_seek(params.shots, rank * params.nt * params.ng * sizeof(float), SEEK_SET);

    /**
     * calculate local objective function & derr & illum & g1(gradient)
     */
    float obj = cal_obj_derr_illum_grad(params, &derr[0], &illum[0], &g1[0], &vv[0], &wlt[0], &sxz[0], &gxz[0]);

    /* MPI reduce */
    MpiInplaceReduce(&obj, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MpiInplaceReduce(&g1[0], params.nz * params.nx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MpiInplaceReduce(&illum[0], params.nz * params.nx, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    objval[iter] = iter == 0 ? obj0 = obj, 1.0 : obj / obj0;

    float epsil = 0;
    float beta = 0;
    if (rank == 0) {
      scale_gradient(&g1[0], &vv[0], &illum[0], params.nz, params.nx, params.precon);
      sf_floatwrite(&illum[0], params.nz * params.nx, params.illums);
      bell_smoothz(&g1[0], &illum[0], params.rbell, params.nz, params.nx);
      bell_smoothx(&illum[0], &g1[0], params.rbell, params.nz, params.nx);
      sf_floatwrite(&g1[0], params.nz * params.nx, params.grads);

      beta = iter == 0 ? 0.0 : cal_beta(&g0[0], &g1[0], &cg[0], params.nz, params.nx);

      cal_conjgrad(&g1[0], &cg[0], beta, params.nz, params.nx);
      epsil = cal_epsilon(&vv[0], &cg[0], params.nz, params.nx);
      cal_vtmp(&vtmp[0], &vv[0], &cg[0], epsil, params.nz, params.nx);

      std::swap(g1, g0); // let g0 be the previous gradient
    }

    sf_seek(params.shots, rank * params.nt * params.ng * sizeof(float), SEEK_SET);
    MPI_Bcast(&vtmp[0], params.nz * params.nx, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float alpha = calVelUpdateStepLen(params, &vtmp[0], &wlt[0], &sxz[0], &gxz[0], &derr[0], epsil);

    if (rank == 0) {
      update_vel(&vv[0], &cg[0], alpha, params.nz, params.nx);

      sf_floatwrite(&vv[0], params.nz * params.nx, params.vupdates);

      // output important information at each FWI iteration
      INFO() << format("iteration %d obj=%f  beta=%f  epsil=%f  alpha=%f") % (iter + 1) % obj % beta % epsil % alpha;
      INFO() << timer.format(2);
    } // end of rank 0

    MPI_Bcast(&vv[0], params.nz * params.nx, MPI_FLOAT, 0, MPI_COMM_WORLD);
  } /// end of iteration

  if (rank == 0) {
    sf_floatwrite(&objval[0], params.niter, params.objs);
  }

  sf_close();
  MPI_Finalize();

  return 0;
}
