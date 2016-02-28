/*
 * spongabc4d.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include <ostream>
#include "spongabc4d.h"
#include "sum.h"

SpongAbc4d::SpongAbc4d(const Velocity &_vel, const FmParams &_params)
  : vel(_vel), params(_params), bndr(params.nb)
{
  initCoeff();
  transVel();
  initbndr();
}

void SpongAbc4d::stepForward(float *p0, float *p1) {
  int nx = vel.nx;
  int nz = vel.nz;
  const std::vector<float> &vv = vel.dat;

//  fprintf(stderr, "before loop, sum p1: %.40f\n", sum(p1, nx * nz));
//  fprintf(stderr, "before loop, p1[3442] %.20f\n", p1[3442]);

  for (int ix = 2; ix < nx - 2; ix++) {
    for (int iz = 2; iz < nz - 2; iz++) {
      float tmp = c0 * p1[ix * nz + iz] +
                  c11 * (p1[ix * nz + (iz - 1)] + p1[ix * nz + (iz + 1)]) +
                  c12 * (p1[ix * nz + (iz - 2)] + p1[ix * nz + (iz + 2)]) +
                  c21 * (p1[(ix - 1) * nz + iz] + p1[(ix + 1) * nz + iz]) +
                  c22 * (p1[(ix - 2) * nz + iz] + p1[(ix + 2) * nz + iz]);

//      if (ix ==41 && iz == 2) {
//        fprintf(stderr, "c0 %.20f\n", c0);
//        fprintf(stderr, "c11 %.20f\n", c11);
//        fprintf(stderr, "c12 %.20f\n", c12);
//        fprintf(stderr, "c21 %.20f\n", c21);
//        fprintf(stderr, "c22 %.20f\n", c22);
//        fprintf(stderr, "p1[ix][iz] %.20f\n", p1[ix * nz + iz]);
//        fprintf(stderr, "p1[ix][iz-1] %.20f\n", p1[ix * nz + (iz - 1)]);
//        fprintf(stderr, "p1[ix][iz+1] %.20f\n", p1[ix * nz + (iz + 1)]);
//        fprintf(stderr, "p1[ix][iz-2] %.20f\n", p1[ix * nz + (iz - 2)]);
//        fprintf(stderr, "p1[ix][iz+2] %.20f\n", p1[ix * nz + (iz + 2)]);
//        fprintf(stderr, "p1[ix-1][iz] %.20f\n", p1[(ix - 1) * nz + iz]);
//        fprintf(stderr, "p1[ix+1][iz] %.20f\n", p1[(ix + 1) * nz + iz]);
//        fprintf(stderr, "p1[ix-2][iz] %.20f\n", p1[(ix - 2) * nz + iz]);
//        fprintf(stderr, "p1[ix+2][iz] %.20f\n", p1[(ix + 2) * nz + iz]);
//        fprintf(stderr, "p1[%d] %.20f\n", (ix + 2) * nz + iz, p1[(ix + 2) * nz + iz]);
//      }
//      fprintf(stderr, "ix %d, iz %d, tmp %.20f\n", ix, iz, tmp);
      p0[ix * nz + iz] = 2 * p1[ix * nz + iz] - p0[ix * nz + iz] + vv[ix * nz + iz] * tmp;
    }
  }

}

void SpongAbc4d::transVel() {
  int nx = vel.nx;
  int nz = vel.nz;
  float dt = params.dt;
  std::vector<float> &vv = vel.dat;

  for(int ix=0;ix<nx;ix++){
    for(int iz=0;iz<nz;iz++){
      float tmp=vv[ix * nz + iz]*dt;
      vv[ix * nz + iz]=tmp*tmp;/* vv=vv^2*dt^2 */
    }
  }
}

void SpongAbc4d::initCoeff() {
  /*< initialize 4-th order FD coefficients >*/
  float dz = params.dz;
  float dx = params.dx;

  float tmp = 1.0/(dz*dz);
  c11 = 4.0*tmp/3.0;
  c12= -tmp/12.0;
  tmp = 1.0/(dx*dx);
  c21 = 4.0*tmp/3.0;
  c22= -tmp/12.0;
  c0=-2.0*(c11+c12+c21+c22);
}

void SpongAbc4d::applySponge(float* p) {
  int nb = params.nb;
  int nz = vel.nz;
  int nx = vel.nx;

  for(int ib=0; ib<nb; ib++) {
    float w = bndr[ib];

    int ibz = nz-ib-1;
    for(int ix=0; ix<nx; ix++) {
      p[ix * nz + ibz] *= w; /* bottom sponge */
    }

    int ibx = nx-ib-1;
    for(int iz=0; iz<nz; iz++) {
      p[ib  * nz + iz] *= w; /*   left sponge */
      p[ibx * nz + iz] *= w; /*  right sponge */
    }
  }
}

void SpongAbc4d::initbndr() {
  int nb = params.nb;

  for(int ib=0;ib<nb;ib++){
    float tmp=0.015*(nb-ib);
    bndr[ib]=expf(-tmp*tmp);
  }

  static bool init = false;
  if (!init) {
    init = true;
    fprintf(stderr, "sum bndr %.20f\n", sum(bndr));
  }
}
