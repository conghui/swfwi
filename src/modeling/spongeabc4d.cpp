/*
 * spongabc4d.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include "spongeabc4d.h"

#include <ostream>
#include <cmath>
#include <functional>
#include "sum.h"
#include <cstdio>
#include "logger.h"

static void expand(Velocity &exvel, const Velocity &v0, int nb) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nxpad = nx + 2 * nb;
  int nzpad = nz + nb;
  const std::vector<float> &a = v0.dat;
  std::vector<float> &b = exvel.dat;

  /// internal
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[(nb + ix) * nzpad + iz] = a[ix * nz + iz];
    }
  }

  /// boundary
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + (nzpad - iz - 1)] = b[ix * nzpad + (nzpad - nb - 1)];  /* bottom*/
    }
  }

  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      b[ix * nzpad + iz] = b[nb * nzpad + iz];/* left */
      b[(nxpad - ix - 1) * nzpad + iz] = b[(nxpad - nb - 1) * nzpad + iz]; /* right */
    }
  }
}


static inline float velTransform(float vv, float dt) {
  float tmp=vv*dt;
  return tmp*tmp;/* vv=vv^2*dt^2 */
}

SpongeAbc4d::SpongeAbc4d(float _dt, float _dx, float _dz, int _nb)
  : IModeling(), bndr(_nb, 0), dt(_dt), dx(_dx), dz(_dz), nb(_nb)
{
  initCoeff();
  initbndr();
}

void SpongeAbc4d::stepForward(float *p0, const float *p1) const {
  int nx = vel->nx;
  int nz = vel->nz;
//  DEBUG() << format("in sponge, nx %d, nz %d") % nx % nz;
  const std::vector<float> &vv = vel->dat;

#pragma omp parallel for
  for (int ix = 2; ix < nx - 2; ix++) {
    for (int iz = 2; iz < nz - 2; iz++) {
      float tmp = c0 * p1[ix * nz + iz] +
                  c11 * (p1[ix * nz + (iz - 1)] + p1[ix * nz + (iz + 1)]) +
                  c12 * (p1[ix * nz + (iz - 2)] + p1[ix * nz + (iz + 2)]) +
                  c21 * (p1[(ix - 1) * nz + iz] + p1[(ix + 1) * nz + iz]) +
                  c22 * (p1[(ix - 2) * nz + iz] + p1[(ix + 2) * nz + iz]);

      float transvel = velTransform(vv[ix * nz + iz], dt);
      p0[ix * nz + iz] = 2 * p1[ix * nz + iz] - p0[ix * nz + iz] + transvel * tmp;
    }
  }

}

void SpongeAbc4d::initCoeff() {
  /*< initialize 4-th order FD coefficients >*/
  float tmp = 1.0/(dz*dz);
  c11 = 4.0*tmp/3.0;
  c12= -tmp/12.0;
  tmp = 1.0/(dx*dx);
  c21 = 4.0*tmp/3.0;
  c22= -tmp/12.0;
  c0=-2.0*(c11+c12+c21+c22);
}

void SpongeAbc4d::applySponge(float* p) const {
  int nz = vel->nz;
  int nx = vel->nx;

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



void SpongeAbc4d::initbndr() {
  for(int ib=0;ib<nb;ib++){
    float tmp=0.015*(nb-ib);
    bndr[ib]=std::exp(-tmp*tmp);
  }
}



Velocity SpongeAbc4d::expandDomain(const Velocity& v0) const {
  int nxpad = v0.nx + 2 * nb;
  int nzpad = v0.nz + nb; // free surface
  Velocity exvel(nxpad, nzpad);

  expand(exvel, v0, nb);
  return exvel;
}


void SpongeAbc4d::addSource(float* p, const float* source, const ShotPosition& pos) const {
  manipSource(p, source, pos, std::plus<float>());
}

void SpongeAbc4d::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) const {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + nb;
    int gz = geoPos.getz(ig);
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }

}

std::vector<float> SpongeAbc4d::initBndryVector(int nt) const {
  if (vel == NULL) {
    ERROR() << __PRETTY_FUNCTION__ << ": you should bind velocity first";
    exit(1);
  }
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nx = nxpad - 2*nb;
  int nz = nzpad - nb;

  bndrSize = FDLEN * (
             nx +   /* bottom */
             2 * nz /* left + right */
             );
  return std::vector<float>(nt*bndrSize, 0);
}

void SpongeAbc4d::writeBndry(float* _bndr, const float* p, int it) const {
  /**
   * say the FDLEN = 2, then the boundary we should save is mark by (*)
   * we omit the upper layer
   *
   *    **+-------------+**
   *    **-             -**
   *    **-             -**
   *    **-             -**
   *    **-             -**
   *    **-             -**
   *    **+-------------+**
   *      ***************
   *      ***************
   *
   */
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nx = nxpad - 2 * nb;
  int nz = nzpad - nb;
  float *bndr = &_bndr[it * bndrSize];

  for (int ix = 0; ix < nx; ix++) {
    for(int iz = 0; iz < FDLEN; iz++) {
      bndr[iz + FDLEN*ix] = p[(ix+nb)*nzpad + (iz+nz)]; // bottom
    }
  }

  for (int iz = 0; iz < nz; iz++) {
    for(int ix=0; ix < FDLEN; ix++) {
      bndr[FDLEN*nx+iz+nz*ix]         = p[(ix-2+nb)*nzpad + (iz)];   // left
      bndr[FDLEN*nx+iz+nz*(ix+FDLEN)] = p[(ix+nx+nb)*nzpad + (iz)];  // right
    }
  }
}

void SpongeAbc4d::readBndry(const float* _bndr, float* p, int it) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nx = nxpad - 2 * nb;
  int nz = nzpad - nb;
  const float *bndr = &_bndr[it * bndrSize];

  for (int ix = 0; ix < nx; ix++) {
    for(int iz = 0; iz < FDLEN; iz++) {
      p[(ix+nb)*nzpad + (iz+nz)] = bndr[iz + FDLEN*ix]; // bottom
    }
  }

  for (int iz = 0; iz < nz; iz++) {
    for(int ix=0; ix < FDLEN; ix++) {
      p[(ix-2+nb)*nzpad + (iz)] = bndr[FDLEN*nx+iz+nz*ix];   // left
      p[(ix+nx+nb)*nzpad + (iz)] = bndr[FDLEN*nx+iz+nz*(ix+FDLEN)];  // right
    }
  }
}

void SpongeAbc4d::stepBackward(float* illum, float* lap, float* p0,
    const float* p1) const {
  int nx = vel->nx;
  int nz = vel->nz;
  const std::vector<float> &vv = vel->dat;

#pragma omp parallel for
  for (int ix = 2; ix < nx - 2; ix++) {
    for (int iz = 2; iz < nz - 2; iz++) {
      float tmp = c0 * p1[ix * nz + iz] +
                  c11 * (p1[ix * nz + (iz - 1)] + p1[ix * nz + (iz + 1)]) +
                  c12 * (p1[ix * nz + (iz - 2)] + p1[ix * nz + (iz + 2)]) +
                  c21 * (p1[(ix - 1) * nz + iz] + p1[(ix + 1) * nz + iz]) +
                  c22 * (p1[(ix - 2) * nz + iz] + p1[(ix + 2) * nz + iz]);

      float transvel = velTransform(vv[ix * nz + iz], dt);
      p0[ix * nz + iz] = 2 * p1[ix * nz + iz] - p0[ix * nz + iz] + transvel * tmp;

      lap[ix * nz + iz] = -4.0 * p1[ix * nz + iz] +
                    p1[ix * nz + (iz - 1)] + p1[ix * nz + (iz + 1)] +
                    p1[(ix - 1) * nz + iz] + p1[(ix + 1) * nz + iz];

//      lap[ix * nz + iz] = -8.0 * p1[ix * nz + iz] +
//                    p1[ix * nz + (iz - 1)] + p1[ix * nz + (iz - 2)] +
//                    p1[ix * nz + (iz + 1)] + p1[ix * nz + (iz + 2)] +
//                    p1[(ix - 1) * nz + iz] + p1[(ix - 2) * nz + iz] +
//                    p1[(ix + 1) * nz + iz] + p1[(ix + 2) * nz + iz];

      illum[ix * nz + iz] += p1[ix * nz + iz] * p1[ix * nz + iz];
    }
  }
}

void SpongeAbc4d::subSource(float* p, const float* source,
    const ShotPosition& pos) const {
  manipSource(p, source, pos, std::minus<float>());
}

void SpongeAbc4d::manipSource(float* p, const float* source,
    const ShotPosition& pos, boost::function2<float, float, float> op) const {
  int nzpad = vel->nz;
  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + nb;
    int sz = pos.getz(is);
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
  }
}

const Velocity& SpongeAbc4d::getVelocity() const {
  return *vel;
}

void SpongeAbc4d::shrinkDomain(float* dst, const float* src, int shrinkNx,
    int shrinkNz) const {
  int nx = shrinkNx;
  int nz = shrinkNz;
  int nzpad = nz + nb;

  /// internal
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      dst[ix * nz + iz] = src[(nb + ix) * nzpad + iz];
    }
  }
}
