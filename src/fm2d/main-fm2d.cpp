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

#include <boost/timer/timer.hpp>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "logger.h"
#include "fm-params.h"
#include "sum.h"
#include "ricker-wavelet.h"
#include "velocity.h"
#include "sf-velocity-reader.h"
#include "common.h"
#include "shot-position.h"

#include "spongabc4d.h"
#include "damp4t10d.h"

int main(int argc, char* argv[])
{
  boost::timer::auto_cpu_timer t;

  /* initialize Madagascar */
  sf_init(argc,argv);

  FmParams &params = FmParams::instance();
  Logger::instance().init("fm");

  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float fm = params.fm;

  std::vector<int> sxz(params.ns); /* source positions */
  std::vector<int> gxz(params.ng); /* geophone positions */

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);
  for (size_t i = 0; i < v0.dat.size(); i++) {
//    DEBUG() << format("%f") % v0.dat[i];
  }
  sf_file f = sf_output("readvel.rsf");
  sf_putint(f, "n1", nz);
  sf_putint(f, "n2", nx);
  sf_floatwrite(&v0.dat[0], nx * nz, f);

  SpongAbc4d fmMethod(dt, params.dx, params.dz, nb);
//  Damp4t10d fmMethod(dt, params.dx, nb);
  Velocity exvel = fmMethod.transformVelocityForModeling(v0);

  fmMethod.setVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  sg_init(&sxz[0], params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  sg_init(&gxz[0], params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);

//  DEBUG() << format("nt %d, fm %d, dt %f, amp %f") % nt % fm % dt % params.amp;
//  DEBUG() << format("sum wlt: %.20f") % sum(wlt);
//
//  DEBUG() << format("sum v0: %.20f") % sum(v0.dat);
//  DEBUG() << format("sum exvel: %.20f") % sum(exvel.dat);
//
//  DEBUG() << "use shot-position";
//  ShotPosition(int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz);
  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);

//  return 0;
  for(int is=0; is<ns; is++)
  {
    std::vector<float> p0(exvel.nz * exvel.nx, 0);
    std::vector<float> p1(exvel.nz * exvel.nx, 0);
    std::vector<float> dobs(params.nt * params.ng, 0);
    std::vector<float> trans(params.nt * params.ng, 0);
    ShotPosition curSrcPos = allSrcPos.clip(is, is);

    for(int it=0; it<nt; it++)
    {
      fmMethod.addSource(&p1[0], &wlt[it], curSrcPos);
      fmMethod.stepForward(&p0[0], &p1[0]);

      std::swap(p1, p0);
//      fprintf(stderr, "is %d, it %d, after swap sum p0: %.20f\n", is , it , sum(p0));
//      fprintf(stderr, "is %d, it %d, after swap sum p1: %.20f\n", is , it , sum(p1));

      record_seis(&dobs[it*ng], &gxz[0], &p0[0], ng, nz, nb);
    }
    matrix_transpose(&dobs[0], &trans[0], ng, nt);
    sf_floatwrite(&trans[0], ng*nt, params.shots);
//    fprintf(stderr, "sum trans: %.20f\n", sum(trans));
  }

  return 0;
}

