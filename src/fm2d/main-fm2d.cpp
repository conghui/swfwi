
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
#include "fm-params.h"
#include "common.h"
#include "ricker-wavelet.h"
#include "cycle-swap.h"
#include "sf-velocity-reader.h"
#include "velocity.h"
#include "spongabc4d.h"
#include "sum.h"

int main(int argc, char **argv) {
  /* initialize Madagascar */
  sf_init(argc, argv);

  Logger::instance().init("fm2d");

  FmParams &params = FmParams::instance();

  // original velocity
  std::vector<float> v0(params.nx * params.nz);
  std::vector<float> wlt(params.nt); /* ricker wavelet */
  std::vector<int> sxz(params.ns); /* source positions */
  std::vector<int> gxz(params.ng); /* geophone positions */

  SfVelocityReader velReader(params.vinit);
  velReader.read(&v0[0], v0.size());
  Velocity vel0(v0, params.nx, params.nz);

  // expand velocity
  Velocity exvel(params.nx + 2 * params.nb, params.nz + params.nb);
  expand(exvel, vel0, params.nb);

  DEBUG() << format("sum_exvel %f") % sum(exvel.dat);

  // create wavelet
  rickerWavelet(&wlt[0], params.nt, params.fm, params.dt);
  sg_init(&sxz[0], params.szbeg, params.sxbeg, params.jsz, params.jsx, params.ns, params.nz);
  sg_init(&gxz[0], params.gzbeg, params.gxbeg, params.jgz, params.jgx, params.ng, params.nz);

  for(int is=0; is < params.ns; is++)
  {
    std::vector<float> p0(exvel.nx * exvel.nz, 0);
    std::vector<float> p1(exvel.nx * exvel.nz, 0);
    std::vector<float> dobs(params.nt * params.ng);
    std::vector<float> trans(params.nt * params.ng);

    SpongAbc4d fm(exvel, params);

      char fn[128];
      sprintf(fn, "fmwld%d.rsf", is);
      sf_file f = sf_output(fn);
      sf_putint(f, "n1", exvel.nz);
      sf_putint(f, "n2", exvel.nx);
      sf_putint(f, "n3", params.nt);

    for(int it=0; it<params.nt; it++)
    {
//      fprintf(stderr, "is %d, it %d, before add source sum p0: %.20f\n", is , it , sum(p0));
//      fprintf(stderr, "is %d, it %d, before add source sum p1: %.20f\n", is , it , sum(p1));

      add_source(&p1[0], &wlt[it], &sxz[is], 1, params.nz, params.nb, true);
//      fprintf(stderr, "is %d, it %d, after add source sum p0: %.20f\n", is , it , sum(p0));
//      fprintf(stderr, "is %d, it %d, after add source sum p1: %.20f\n", is , it , sum(p1));

      fm.stepForward(&p0[0], &p1[0]);
//      fprintf(stderr, "is %d, it %d, after forward sum p0: %.20f\n", is , it , sum(p0));
//      fprintf(stderr, "is %d, it %d, after forward sum p1: %.20f\n", is , it , sum(p1));

      sf_floatwrite(&p0[0], p0.size(), f);

      fm.applySponge(&p0[0]);
      fm.applySponge(&p1[0]);

      std::swap(p0, p1);
//      fprintf(stderr, "is %d, it %d, after swap sum p0: %.20f\n", is , it , sum(p0));
//      fprintf(stderr, "is %d, it %d, after swap sum p1: %.20f\n", is , it , sum(p1));

      record_seis(&dobs[it * params.ng], &gxz[0], &p0[0], params.ng, params.nz, params.nb);
    }
    matrix_transpose(&dobs[0], &trans[0], params.ng, params.nt);
    sf_floatwrite(&trans[0], params.ng*params.nt, params.shots);
  }

  return 0;
}
