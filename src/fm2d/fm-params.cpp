/*
 * fm-params.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include "fm-params.h"


#include <cmath>

FmParams *FmParams::ins = NULL;

FmParams::FmParams() {
  /*< set up I/O files >*/
  vinit=sf_input ("vinit");   /* initial velocity model, unit=m/s */
  shots=sf_output("shots");
}

FmParams::~FmParams() {
  sf_close();
  delete ins;
}

FmParams& FmParams::instance() {
  if (ins == NULL) {
    ins = new FmParams();
    ins->getInputParams();
    ins->check();
    ins->putOutputParams();
  }

  return *ins;
}

void FmParams::getInputParams() {
  /* get parameters for forward modeling */
    if (!sf_histint(vinit,"n1",&nz)) sf_error("no n1");
    if (!sf_histint(vinit,"n2",&nx)) sf_error("no n2");
    if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
    if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

    if (!sf_getfloat("amp",&amp)) amp=1000;
    /* maximum amplitude of ricker */
    if (!sf_getfloat("fm",&fm)) fm=10;
    /* dominant freq of ricker */
    if (!sf_getint("nb",&nb))   nb=30;
    /* thickness of sponge ABC  */
    if (!sf_getfloat("dt",&dt)) sf_error("no dt");
    /* time interval */
    if (!sf_getint("nt",&nt))   sf_error("no nt");
    /* total modeling time steps */
    if (!sf_getint("ns",&ns))   sf_error("no ns");
    /* total shots */
    if (!sf_getint("ng",&ng))   sf_error("no ng");
    /* total receivers in each shot */
    if (!sf_getint("jsx",&jsx))   sf_error("no jsx");
    /* source x-axis  jump interval  */
    if (!sf_getint("jsz",&jsz))   jsz=0;
    /* source z-axis jump interval  */
    if (!sf_getint("jgx",&jgx))   jgx=1;
    /* receiver x-axis jump interval */
    if (!sf_getint("jgz",&jgz))   jgz=0;
    /* receiver z-axis jump interval */
    if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");
    /* x-begining index of sources, starting from 0 */
    if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");
    /* z-begining index of sources, starting from 0 */
    if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");
    /* x-begining index of receivers, starting from 0 */
    if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");
    /* z-begining index of receivers, starting from 0 */
}

void FmParams::putOutputParams() {
  sf_putint(shots,"n1",nt);
  sf_putint(shots,"n2",ng);
  sf_putint(shots,"n3",ns);
  sf_putfloat(shots,"d1",dt);
  sf_putfloat(shots,"d2",jgx*dx);
  sf_putfloat(shots,"o1",0);
  sf_putfloat(shots,"o2",0);
  sf_putstring(shots,"label1","Time");
  sf_putstring(shots,"label2","Lateral");
  sf_putstring(shots,"label3","Shot");
  sf_putstring(shots,"unit1","sec");
  sf_putstring(shots,"unit2","m");
  sf_putfloat(shots,"amp",amp);
  sf_putfloat(shots,"fm",fm);
  sf_putint(shots,"ng",ng);
  sf_putint(shots,"szbeg",szbeg);
  sf_putint(shots,"sxbeg",sxbeg);
  sf_putint(shots,"gzbeg",gzbeg);
  sf_putint(shots,"gxbeg",gxbeg);
  sf_putint(shots,"jsx",jsx);
  sf_putint(shots,"jsz",jsz);
  sf_putint(shots,"jgx",jgx);
  sf_putint(shots,"jgz",jgz);
}

void FmParams::check() {
  if (!(sxbeg >= 0 && szbeg >= 0 && sxbeg + (ns - 1)*jsx < nx && szbeg + (ns - 1)*jsz < nz)) {
    sf_warning("sources exceeds the computing zone!\n");
    exit(1);
  }

  if (!(gxbeg >= 0 && gzbeg >= 0 && gxbeg + (ng - 1)*jgx < nx && gzbeg + (ng - 1)*jgz < nz)) {
    sf_warning("geophones exceeds the computing zone!\n");
    exit(1);
  }

}
