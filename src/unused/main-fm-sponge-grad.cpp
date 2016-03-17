
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

#include "zjh4t10dsponge.h"
#include "sfutil.h"

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

  ShotPosition allSrcPos(params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  ShotPosition allGeoPos(params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);
  Zjh4t10dSponge fmMethod(allSrcPos, allGeoPos, dt, params.dx, params.fm, nb, nt);

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);
  sfFloatWrite2d("v0.rsf", &v0.dat[0], nz, nx);


  Velocity exvel = fmMethod.expandDomain(v0);
  sfFloatWrite2d("exvel.rsf", &exvel.dat[0], exvel.nz, exvel.nx);

  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);


  for(int is=0; is<ns; is++) {
    boost::timer::cpu_timer timer;

    std::vector<float> p0(exvel.nz * exvel.nx, 0);
    std::vector<float> p1(exvel.nz * exvel.nx, 0);
    std::vector<float> dobs(params.nt * params.ng, 0);
    std::vector<float> trans(params.nt * params.ng, 0);
    ShotPosition curSrcPos = allSrcPos.clipRange(is, is);

    for(int it=0; it<nt; it++) {
      fmMethod.addSource(&p1[0], &wlt[it], curSrcPos);
      fmMethod.stepForward(&p0[0], &p1[0]);

      std::swap(p1, p0);
      fmMethod.recordSeis(&dobs[it*ng], &p0[0]);
    }

    matrix_transpose(&dobs[0], &trans[0], ng, nt);
    sf_floatwrite(&trans[0], ng*nt, params.shots);

    DEBUG() << format("shot %d, time %s") % is % timer.format(2).c_str();
  }

  return 0;
}

