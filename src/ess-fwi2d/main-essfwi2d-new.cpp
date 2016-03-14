extern "C" {
#include <rsf.h>
}

#include <cstdlib>

#include "logger.h"
#include "essfwi-params.h"
#include "shot-position.h"
#include "damp4t10d.h"
#include "sf-velocity-reader.h"
#include "ricker-wavelet.h"
#include "essfwiframework.h"
#include "shotdata-reader.h"
#include "sum.h"
#include "updatevelop.h"

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
  Damp4t10d fmMethod(allSrcPos, allGeoPos, dt, dx, fm, nb, nt);

  SfVelocityReader velReader(params.vinit);
  Velocity v0 = SfVelocityReader::read(params.vinit, nx, nz);
  Velocity exvel = fmMethod.expandDomain(v0);
  fmMethod.bindVelocity(exvel);

  std::vector<float> wlt(nt);
  rickerWavelet(&wlt[0], nt, fm, dt, params.amp);

  std::vector<float> dobs(ns * nt * ng);     /* all observed data */
  ShotDataReader::serialRead(params.shots, &dobs[0], ns, nt, ng);

  float vmax = 5500;
  float vmin = 1500;
  UpdateVelOp updatevelop(vmin, vmax, dx, dt);

  int max_iter_update_alpha = 5;
  float maxdv = 200;
  UpdateSteplenOp updateSteplenOp(fmMethod, updatevelop, max_iter_update_alpha, maxdv);

  EssFwiFramework essfwi(fmMethod, updateSteplenOp, updatevelop, wlt, dobs);
  for (int iter = 0; iter < params.niter; iter++) {
    essfwi.epoch(iter);
    essfwi.writeVel(params.vupdates);
  } /// end of iteration

  sf_close();

  return 0;
}
