extern "C" {
#include <rsf.h>
}

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cstdlib>
#include <cmath>

#include "logger.h"
#include "shot-position.h"
#include "damp4t10d.h"
#include "sf-velocity-reader.h"
#include "ricker-wavelet.h"
#include "fwiframework.h"
#include "shotdata-reader.h"
#include "updatevelop.h"
#include "environment.h"

namespace {
class Params {
public:
  Params();
  ~Params();
  void check();

public:
  sf_file vinit;        /* initial velocity model, unit=m/s */
  sf_file shots;        /* recorded shots from exact velocity model */
  sf_file vupdates;     /* updated velocity in iterations */
  sf_file absobjs;         /* absolute values of objective function in iterations */
  sf_file norobjs;         /* normalize values of objective function in iterations */
  int niter;            /* # of iterations */
  int nb;               /* size of the boundary */
  float vmin;
  float vmax;
  float maxdv;
  int nita;
  int seed;

public: // parameters from input files
  int nz;
  int nx;
  float dz;
  float dx;
  int nt;
  int ng;
  int ns;
  float dt;
  float amp;
  float fm;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  int jsx;
  int jsz;
  int jgx;
  int jgz;
};

Params::Params() {
  vinit = sf_input ("vin");       /* initial velocity model, unit=m/s */
  shots = sf_input("shots");      /* recorded shots from exact velocity model */
  vupdates = sf_output("vout");   /* updated velocity in iterations */
  absobjs = sf_output("absobjs"); /* absolute values of objective function in iterations */
  norobjs = sf_output("norobjs"); /* normalized values of objective function in iterations */

  if (!sf_getint("niter", &niter)) { sf_error("no niter"); }      /* number of iterations */
  if (!sf_getfloat("maxdv", &maxdv)) sf_error("no maxdv");        /* max delta v update two iteration*/
  if (!sf_getint("nita", &nita))   { sf_error("no nita"); }       /* max iter refining alpha */
  if (!sf_getint("seed", &seed))   { seed = 10; }                 /* seed for random numbers */

  /* get parameters from velocity model and recorded shots */
  if (!sf_histint(vinit, "n1", &nz)) { sf_error("no n1"); }       /* nz */
  if (!sf_histint(vinit, "n2", &nx)) { sf_error("no n2"); }       /* nx */
  if (!sf_histfloat(vinit, "d1", &dz)) { sf_error("no d1"); }     /* dz */
  if (!sf_histfloat(vinit, "d2", &dx)) { sf_error("no d2"); }     /* dx */
  if (!sf_histint(shots, "n1", &nt)) { sf_error("no nt"); }       /* total modeling time steps */
  if (!sf_histint(shots, "n2", &ng)) { sf_error("no ng"); }       /* total receivers in each shot */
  if (!sf_histint(shots, "n3", &ns)) { sf_error("no ns"); }       /* number of shots */
  if (!sf_histfloat(shots, "d1", &dt)) { sf_error("no dt"); }     /* time sampling interval */
  if (!sf_histfloat(shots, "amp", &amp)) { sf_error("no amp"); }  /* maximum amplitude of ricker */
  if (!sf_histfloat(shots, "fm", &fm)) { sf_error("no fm"); }     /* dominant freq of ricker */
  if (!sf_histint(shots, "sxbeg", &sxbeg)) {sf_error("no sxbeg");} /* x-begining index of sources, starting from 0 */
  if (!sf_histint(shots, "szbeg", &szbeg)) {sf_error("no szbeg");} /* x-begining index of sources, starting from 0 */
  if (!sf_histint(shots, "gxbeg", &gxbeg)) {sf_error("no gxbeg");} /* x-begining index of receivers, starting from 0 */
  if (!sf_histint(shots, "gzbeg", &gzbeg)) {sf_error("no gzbeg");} /* x-begining index of receivers, starting from 0 */
  if (!sf_histint(shots, "jsx", &jsx)) {sf_error("no jsx"); }       /* source x-axis  jump interval  */
  if (!sf_histint(shots, "jsz", &jsz)) { sf_error("no jsz"); }      /* source z-axis jump interval  */
  if (!sf_histint(shots, "jgx", &jgx)) { sf_error("no jgx"); }      /* receiver x-axis jump interval  */
  if (!sf_histint(shots, "jgz", &jgz)) { sf_error("no jgz"); }      /* receiver z-axis jump interval  */
  if (!sf_histint(shots,  "nb",&nb))        { sf_error("no nb"); }  /* thickness of sponge ABC  */
  if (!sf_histfloat(shots, "vmin", &vmin)) { sf_error("no vmin"); } /* minimal velocity in real model*/
  if (!sf_histfloat(shots, "vmax", &vmax)) { sf_error("no vmax"); } /* maximal velocity in real model*/


  /**
   * output parameters
   */
  sf_putint(vupdates, "n1", nz);
  sf_putint(vupdates, "n2", nx);
  sf_putfloat(vupdates, "d1", dz);
  sf_putfloat(vupdates, "d2", dx);
  sf_putstring(vupdates, "label1", "Depth");
  sf_putstring(vupdates, "label2", "Distance");
  sf_putstring(vupdates, "label3", "Iteration");
  sf_putint(vupdates, "n3", niter);
  sf_putint(vupdates, "d3", 1);
  sf_putint(vupdates, "o3", 1);
  sf_putint(absobjs, "n1", niter + 1);
  sf_putfloat(absobjs, "d1", 1);
  sf_putfloat(absobjs, "o1", 1);
  sf_putstring(absobjs, "label1", "Absolute");
  sf_putint(norobjs, "n1", niter + 1);
  sf_putfloat(norobjs, "d1", 1);
  sf_putfloat(norobjs, "o1", 1);
  sf_putstring(norobjs, "label1", "Normalize");

  check();
}

Params::~Params() {
  sf_close();
}

void Params::check() {
  if (!(sxbeg >= 0 && szbeg >= 0 && sxbeg + (ns - 1)*jsx < nx && szbeg + (ns - 1)*jsz < nz)) {
    sf_warning("sources exceeds the computing zone!\n");
    exit(1);
  }

  if (!(gxbeg >= 0 && gzbeg >= 0 && gxbeg + (ng - 1)*jgx < nx && gzbeg + (ng - 1)*jgz < nz)) {
    sf_warning("geophones exceeds the computing zone!\n");
    exit(1);
  }
}

} /// end of name space


int main(int argc, char *argv[]) {
  sf_init(argc, argv);                /* initialize Madagascar */
  Environment::setDatapath();
  Params params;

  /// configure logger
  const char *logfile = "fwi-damp.log";
  FILELog::setLogFile(logfile);

  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float fm = params.fm;
  float dx = params.dx;
  float vmin = params.vmin;
  float vmax = params.vmax;
  float maxdv = params.maxdv;
  int nita = params.nita;

  srand(params.seed);

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

  FwiUpdateVelOp updatevelop(vmin, vmax, dx, dt);
  FwiUpdateSteplenOp updateSteplenOp(fmMethod, updatevelop, nita, maxdv);

  FwiFramework fwi(fmMethod, updateSteplenOp, updatevelop, wlt, dobs);

  std::vector<float> absobj;
  std::vector<float> norobj;
  for (int iter = 0; iter < params.niter; iter++) {
		fwi.epoch(iter);
    fwi.writeVel(params.vupdates);
    float obj = fwi.getUpdateObj();
    if (iter == 0) {
      float obj0 = fwi.getInitObj();
      absobj.push_back(obj0);
      norobj.push_back(1);
    }
    absobj.push_back(obj);
    norobj.push_back(obj / absobj[0]);
  } /// end of iteration

  sf_floatwrite(&absobj[0], absobj.size(), params.absobjs);
  sf_floatwrite(&norobj[0], norobj.size(), params.norobjs);

  sf_close();

  return 0;
}
