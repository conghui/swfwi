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


void expand2d(float** b, float** a, int nx, int nz, int nxpad, int nzpad, int nb)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
  int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz)			\
    shared(b,a,nb,nz,nx)
#endif
  for     (ix=0;ix<nx;ix++) {
    for (iz=0;iz<nz;iz++) {
      b[nb+ix][iz] = a[ix][iz];
    }
  }

  for     (ix=0; ix<nxpad; ix++)
    for (iz=0; iz<nb;    iz++)
      b[ix][nzpad-iz-1] = b[ix][nzpad-nb-1];/* bottom*/

  for     (ix=0; ix<nb;    ix++) {
    for (iz=0; iz<nzpad; iz++) {
      b[ix 	 ][iz] = b[nb  		][iz];/* left */
      b[nxpad-ix-1 ][iz] = b[nxpad-nb-1	][iz];/* right */
    }
  }
}


void window2d(float **a, float **b, int nx, int nz, int nb)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
  int iz,ix;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz)			\
    shared(b,a,nb,nz,nx)
#endif
  for     (ix=0;ix<nx;ix++) {
    for (iz=0;iz<nz;iz++) {
      a[ix][iz]=b[nb+ix][iz] ;
    }
  }
}


void step_forward(float **p0, float **p1, float **vv, float dz, float dx, int nxpad, int nzpad)
/*< forward modeling step >*/
{
  int ix,iz;
  float tmp;
  float c0, c11, c12, c21, c22;

  /**
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
    private(ix,iz,tmp)			\
    shared(p1,p0,vv,nzpad,nxpad,c0,c11,c12,c21,c22)
#endif
*/

  /*< initialize 4-th order FD coefficients >*/
  tmp = 1.0/(dz*dz);
  c11 = 4.0*tmp/3.0;
  c12= -tmp/12.0;
  tmp = 1.0/(dx*dx);
  c21 = 4.0*tmp/3.0;
  c22= -tmp/12.0;
  c0=-2.0*(c11+c12+c21+c22);

  for (ix=2; ix < nxpad-2; ix++)
    for (iz=2; iz < nzpad-2; iz++)
    {
      tmp =	c0*p1[ix][iz]+
          c11*(p1[ix][iz-1]+p1[ix][iz+1])+
          c12*(p1[ix][iz-2]+p1[ix][iz+2])+
          c21*(p1[ix-1][iz]+p1[ix+1][iz])+
          c22*(p1[ix-2][iz]+p1[ix+2][iz]);

      p0[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+vv[ix][iz]*tmp;
    }
}

void apply_sponge(float**p0, const float *bndr, int nxpad, int nzpad, int nb)
/*< apply absorbing boundary condition >*/
{
  int ix,iz,ib,ibx,ibz;
  float w;

#ifdef _OPENMP
#pragma omp parallel for		\
    private(ib,iz,ix,ibz,ibx,w)		\
    shared(p0,bndr,nzpad,nxpad,nb)
#endif
  for(ib=0; ib<nb; ib++) {
    w = bndr[ib];

    ibz = nzpad-ib-1;
    for(ix=0; ix<nxpad; ix++) {
      p0[ix][ibz] *= w; /* bottom sponge */
    }

    ibx = nxpad-ib-1;
    for(iz=0; iz<nzpad; iz++) {
      p0[ib ][iz] *= w; /*   left sponge */
      p0[ibx][iz] *= w; /*  right sponge */
    }
  }
}


void add_source(int *sxz, float **p, int ns, float *source, int nz, int nb, bool add)
/*< add source term >*/
{
  int is, sx, sz;
  if(add){
    for(is=0;is<ns; is++){
      sx=sxz[is]/nz+nb;
      sz=sxz[is]%nz;
      p[sx][sz]+=source[is];
    }
  }else{
    for(is=0;is<ns; is++){
      sx=sxz[is]/nz+nb;
      sz=sxz[is]%nz;
      p[sx][sz]-=source[is];
    }
  }
}

void record_seis(float *seis_it, int *gxz, float **p, int ng, int nz, int nb)
/*< record seismogram at time it into a vector length of ng >*/
{
  int ig, gx, gz;
  for(ig=0;ig<ng; ig++)
  {
    gx=gxz[ig]/nz+nb;
    gz=gxz[ig]%nz;
    seis_it[ig]=p[gx][gz];
  }
}

void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
  int i1, i2;

  for(i2=0; i2<n2; i2++)
    for(i1=0; i1<n1; i1++)
      trans[i2+n2*i1]=matrix[i1+n1*i2];
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz)
/*< shot/geophone position initialize >*/
{
  int is, sz, sx;
  for(is=0; is<ns; is++)
  {
    sz=szbeg+is*jsz;
    sx=sxbeg+is*jsx;
    sxz[is]=sz+nz*sx;
  }
}

int main(int argc, char* argv[])
{
  /* initialize Madagascar */
  sf_init(argc,argv);
#ifdef _OPENMP
  omp_init();
#endif

  const FmParams &params = FmParams::instance();
  int nz = params.nz;
  int nx = params.nx;
  int nb = params.nb;
  int ng = params.ng;
  int nt = params.nt;
  int ns = params.ns;
  float dt = params.dt;
  float dz = params.dz;
  float dx = params.dx;
  float fm = params.fm;
  int nzpad=nz+nb;
  int nxpad=nx+2*nb;

  float **v0=sf_floatalloc2(nz,nx);
  float **vv=sf_floatalloc2(nzpad, nxpad);
  float **p0=sf_floatalloc2(nzpad, nxpad);
  float **p1=sf_floatalloc2(nzpad, nxpad);
  float *dobs=(float*)malloc(ng*nt*sizeof(float));
  float *trans=(float*)malloc(ng*nt*sizeof(float));
  float *bndr=sf_floatalloc(nb);
  int *sxz=sf_intalloc(ns);
  int *gxz=sf_intalloc(ng);

  sf_floatread(v0[0],nz*nx,params.vinit);
  expand2d(vv, v0, nx, nz, nxpad, nzpad, nb);


  memset(dobs,0,ng*nt*sizeof(float));
  memset(trans,0,ng*nt*sizeof(float));
  for(int ib=0;ib<nb;ib++){
    float tmp=0.015*(nb-ib);
    bndr[ib]=expf(-tmp*tmp);
  }
  for(int ix=0;ix<nxpad;ix++){
    for(int iz=0;iz<nzpad;iz++){
      float tmp=vv[ix][iz]*dt;
      vv[ix][iz]=tmp*tmp;/* vv=vv^2*dt^2 */
    }
  }

  std::vector<float> wlt(nt);
//  rickerWavelet(&wlt[0], nt, fm, dt);
  for(int it=0; it<nt; it++){
    float tmp=SF_PI*fm*(it*dt-1.0/fm);tmp=tmp*tmp;
    wlt[it]=params.amp*(1.0-2.0*tmp)*expf(-tmp);
  }

  sg_init(sxz, params.szbeg, params.sxbeg, params.jsz, params.jsx, ns, nz);
  sg_init(gxz, params.gzbeg, params.gxbeg, params.jgz, params.jgx, ng, nz);

  fprintf(stderr, "sum bndr %.20f\n", sum(bndr, nb));
  for(int is=0; is<ns; is++)
  {
    memset(p0[0],0,nzpad*nxpad*sizeof(float));
    memset(p1[0],0,nzpad*nxpad*sizeof(float));

      char fn[128];
      sprintf(fn, "mswld%d.rsf", is);
      sf_file f = sf_output(fn);
      sf_putint(f, "n1", nzpad);
      sf_putint(f, "n2", nxpad);
      sf_putint(f, "n3", nt);

    for(int it=0; it<nt; it++)
    {
      int size = nzpad * nxpad;
      add_source(&sxz[is], p1, 1, &wlt[it], nz, nb, true);
      step_forward(p0, p1, vv, dz, dx, nxpad, nzpad);

      sf_floatwrite(p0[0], size, f);

      apply_sponge(p1, bndr, nxpad, nzpad, nb);
      apply_sponge(p0, bndr, nxpad, nzpad, nb);

      float **ptr=p0; p0=p1; p1=ptr;

      record_seis(&dobs[it*ng], gxz, p0, ng, nz, nb);
    }
    matrix_transpose(dobs, trans, ng, nt);
    sf_floatwrite(trans, ng*nt, params.shots);

  }

  free(sxz);
  free(gxz);
  free(dobs);
  free(trans);
  free(*v0); free(v0);
  free(*vv); free(vv);
  free(*p0); free(p0);
  free(*p1); free(p1);
  free(bndr);

  exit(0);
}

