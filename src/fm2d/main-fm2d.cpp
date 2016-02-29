/* 2-D prestack forward modeling using sponge ABC using 4-th order FD
NB: prepare high quality prestack seismic data for LSRTM and FWI
Top boundary is free surface (no ABC applied)!
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
 */

extern "C" {
#include <rsf.h>
}
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fm-params.h"
#include "sum.h"
#include "ricker-wavelet.h"
#include "velocity.h"
#include "spongabc4d.h"
#include "sf-velocity-reader.h"
#include "common.h"

int main(int argc, char* argv[])
{
  /* initialize Madagascar */
  sf_init(argc,argv);

  FmParams &params = FmParams::instance();
  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float fm = params.fm;

  std::vector<float> v0(nx * nz);
  std::vector<int> sxz(params.ns); /* source positions */
  std::vector<int> gxz(params.ng); /* geophone positions */

  SfVelocityReader velReader(params.vinit);
  velReader.read(&v0[0], v0.size());

  Velocity vv0(v0, nx, nz);
  Velocity exvel(nx+2*nb, nz+nb);
  expand(exvel, vv0, nb);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  sg_init(&sxz[0], params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  sg_init(&gxz[0], params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);

  SpongAbc4d fmMethod(exvel, params);
  for(int is=0; is<ns; is++)
  {
    std::vector<float> p0(exvel.nz * exvel.nx, 0);
    std::vector<float> p1(exvel.nz * exvel.nx, 0);
    std::vector<float> dobs(params.nt * params.ng, 0);
    std::vector<float> trans(params.nt * params.ng, 0);

    for(int it=0; it<nt; it++)
    {
      add_source(&p1[0], &wlt[it], &sxz[is], 1, nz, nb, true);
      fmMethod.stepForward(&p0[0], &p1[0]);

      std::swap(p1, p0);

      record_seis(&dobs[it*ng], &gxz[0], &p0[0], ng, nz, nb);
    }
    matrix_transpose(&dobs[0], &trans[0], ng, nt);
    sf_floatwrite(&trans[0], ng*nt, params.shots);

  }

  return 0;
}

