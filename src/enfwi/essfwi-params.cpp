/*
 * essfwi-params.cpp
 *
 *  Created on: Feb 27, 2016
 *      Author: rice
 */

#include "essfwi-params.h"


#include <cmath>

EssFwiParams *EssFwiParams::ins = NULL;

EssFwiParams::EssFwiParams() {
  vinit = sf_input ("vin");       /* initial velocity model, unit=m/s */
  shots = sf_input("shots");      /* recorded shots from exact velocity model */

  vupdates = sf_output("vout");   /* updated velocity in iterations */
  grads = sf_output("grads");     /* gradient in iterations */
  illums = sf_output("illums");   /* source illumination in iterations */
  objs = sf_output("objs");       /* values of objective function in iterations */
}

EssFwiParams::~EssFwiParams() {
  sf_close();
  delete ins;
}

EssFwiParams& EssFwiParams::instance() {
  if (ins == NULL) {
    ins = new EssFwiParams();
    ins->getInputParams();
    ins->check();
    ins->putOutputParams();
  }

  return *ins;
}

void EssFwiParams::getInputParams() {
  /* get parameters from velocity model and recorded shots */
   if (!sf_getbool("verb", &verb)) {
     verb = true;  /* vebosity */
   }
   if (!sf_histint(vinit, "n1", &nz)) {
     sf_error("no n1");  /* nz */
   }
   if (!sf_histint(vinit, "n2", &nx)) {
     sf_error("no n2");  /* nx */
   }
   if (!sf_histfloat(vinit, "d1", &dz)) {
     sf_error("no d1");  /* dz */
   }
   if (!sf_histfloat(vinit, "d2", &dx)) {
     sf_error("no d2");  /* dx */
   }
   if (!sf_getbool("precon", &precon)) {
     precon = false;  /* precondition or not */
   }
   if (!sf_getint("niter", &niter)) {
     niter = 100;  /* number of iterations */
   }
   if (!sf_getint("rbell", &rbell)) {
     rbell = 2;  /* radius of bell smooth */
   }

   if (!sf_histint(shots, "n1", &nt)) {
     sf_error("no nt");
   }
   /* total modeling time steps */
   if (!sf_histint(shots, "n2", &ng)) {
     sf_error("no ng");
   }
   /* total receivers in each shot */
   if (!sf_histint(shots, "n3", &ns)) {
     sf_error("no ns");
   }
   /* number of shots */
   if (!sf_histfloat(shots, "d1", &dt)) {
     sf_error("no dt");
   }
   /* time sampling interval */
   if (!sf_histfloat(shots, "amp", &amp)) {
     sf_error("no amp");
   }
   /* maximum amplitude of ricker */
   if (!sf_histfloat(shots, "fm", &fm)) {
     sf_error("no fm");
   }
   /* dominant freq of ricker */
   if (!sf_histint(shots, "sxbeg", &sxbeg)) {
     sf_error("no sxbeg");
   }
   /* x-begining index of sources, starting from 0 */
   if (!sf_histint(shots, "szbeg", &szbeg)) {
     sf_error("no szbeg");
   }
   /* x-begining index of sources, starting from 0 */
   if (!sf_histint(shots, "gxbeg", &gxbeg)) {
     sf_error("no gxbeg");
   }
   /* x-begining index of receivers, starting from 0 */
   if (!sf_histint(shots, "gzbeg", &gzbeg)) {
     sf_error("no gzbeg");
   }
   /* x-begining index of receivers, starting from 0 */
   if (!sf_histint(shots, "jsx", &jsx)) {
     sf_error("no jsx");
   }
   /* source x-axis  jump interval  */
   if (!sf_histint(shots, "jsz", &jsz)) {
     sf_error("no jsz");
   }
   /* source z-axis jump interval  */
   if (!sf_histint(shots, "jgx", &jgx)) {
     sf_error("no jgx");
   }
   /* receiver x-axis jump interval  */
   if (!sf_histint(shots, "jgz", &jgz)) {
     sf_error("no jgz");
   }
   /* receiver z-axis jump interval  */

   if (!sf_getint("nb",&nb))   nb=30;
   /* thickness of sponge ABC  */


   /* filename of the shot data */
   if (!(obsDataFileName = sf_histstring(shots, "in"))) {
     sf_error("cannot find observed data file path");
   }

}

void EssFwiParams::putOutputParams() {
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
  sf_putint(grads, "n1", nz);
  sf_putint(grads, "n2", nx);
  sf_putint(grads, "n3", niter);
  sf_putfloat(grads, "d1", dz);
  sf_putfloat(grads, "d2", dx);
  sf_putint(grads, "d3", 1);
  sf_putint(grads, "o3", 1);
  sf_putstring(grads, "label1", "Depth");
  sf_putstring(grads, "label2", "Distance");
  sf_putstring(grads, "label3", "Iteration");
  sf_putint(illums, "n1", nz);
  sf_putint(illums, "n2", nx);
  sf_putfloat(illums, "d1", dz);
  sf_putfloat(illums, "d2", dx);
  sf_putint(illums, "n3", niter);
  sf_putint(illums, "d3", 1);
  sf_putint(illums, "o3", 1);
  sf_putint(objs, "n1", niter);
  sf_putint(objs, "n2", 1);
  sf_putfloat(objs, "d1", 1);
  sf_putfloat(objs, "o1", 1);
}

void EssFwiParams::check() {
  if (!(sxbeg >= 0 && szbeg >= 0 && sxbeg + (ns - 1)*jsx < nx && szbeg + (ns - 1)*jsz < nz)) {
    sf_warning("sources exceeds the computing zone!\n");
    exit(1);
  }

  if (!(gxbeg >= 0 && gzbeg >= 0 && gxbeg + (ng - 1)*jgx < nx && gzbeg + (ng - 1)*jgz < nz)) {
    sf_warning("geophones exceeds the computing zone!\n");
    exit(1);
  }

}
