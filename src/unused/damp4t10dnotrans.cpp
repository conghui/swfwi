/*
 * damp4t10dnotrans.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#include "damp4t10dnotrans.h"


#include <cmath>
#include <functional>
#include "logger.h"
#include "sum.h"
#include "fd4t10s-damp-zjh.h"
#include "fd4t10s-zjh.h"
#include "sfutil.h"
#include "common.h"

extern "C" {
#include <rsf.h>
}

static void initbndr(std::vector<float> &bndr, int nb) {
  for(int ib=0;ib<nb;ib++){
    float tmp=0.015*(nb-ib);
    bndr[ib]=std::exp(-tmp*tmp);
  }
}


static void fillForStencil(Velocity &exvel, int halo) {
  int nxpad = exvel.nx;
  int nzpad = exvel.nz;
  int nx = nxpad - 2 * halo;
  int nz = nzpad - 2 * halo;

  std::vector<float> &vel_e = exvel.dat;

  //expand z direction first
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = 0; iz < halo; iz++) {
      vel_e[ix * nzpad + iz] = vel_e[ix * nzpad + halo];               // top
    }
    for (int iz = nz + halo; iz < nzpad; iz ++) {
      vel_e[ix * nzpad + iz] = vel_e[ix * nzpad + (nz + halo - 1)]; // bottom
    }
  }

  //Then x direction
  for (int iz = 0; iz < nzpad; iz++) {
    for (int ix = 0; ix < halo; ix++) {
      vel_e[ix * nzpad + iz] = vel_e[halo * nzpad + iz];               // left
    }
    for (int ix = nx + halo; ix < nxpad; ix++) {
      vel_e[ix * nzpad + iz] = vel_e[(nx + halo - 1) * nzpad + iz]; // right
    }
  }

}

static void expandForStencil(Velocity &exvel, const Velocity &v0, int halo) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nzpad = nz + 2 * halo;

  const std::vector<float> &vel = v0.dat;
  std::vector<float> &vel_e = exvel.dat;

  //copy the vel into vel_e
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = halo; iz < nz + halo; iz++) {
      vel_e[ix * nzpad + iz] = vel[(ix - halo) * nz +  (iz - halo)];
    }
  }

  fillForStencil(exvel, halo);

}

static void expandBndry(Velocity &exvel, const Velocity &v0, int nb) {
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

  /// boundary, free surface
  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[ix * nzpad + iz] = a[iz];                              /* left */
      b[(nb + nx + ix) * nzpad + iz] = a[(nx - 1) * nz + iz];  /* right */
    }
  }

  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + (nz + iz)] = b[ix * nzpad + (nz - 1)];  /* bottom*/
    }
  }
}


Velocity Damp4t10dNotrans::expandDomain(const Velocity& _vel) {
  // expand for boundary, free surface
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + nb);
  expandBndry(exvelForBndry, _vel, nb);

  // expand for stencil
  Velocity ret(exvelForBndry.nx+2*EXFDBNDRYLEN, exvelForBndry.nz+2*EXFDBNDRYLEN);
  expandForStencil(ret, exvelForBndry, EXFDBNDRYLEN);

  return ret;
}

void Damp4t10dNotrans::stepForward(float* p0, float* p1) const {
  fd4t10s_damp_zjh_2d_vnotrans(p0, p1, &vel->dat[0], vel->nx, vel->nz, nb + EXFDBNDRYLEN, dx, dt);
//  fd4t10s_zjh_2d(p0, p1, &vel->dat[0], vel->nx, vel->nz, nb + FDLEN, dx, dt);
//  applySponge(p0, &bndr[0], vel->nx, vel->nz, nb + FDLEN);
//  applySponge(p1, &bndr[0], vel->nx, vel->nz, nb + FDLEN);
}

void Damp4t10dNotrans::bindVelocity(const Velocity& _vel) {
  this->vel = &_vel;
}

void Damp4t10dNotrans::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) const {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + nb + EXFDBNDRYLEN;
    int gz = geoPos.getz(ig) + EXFDBNDRYLEN;
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }
//  DEBUG() << format("sum %.20f") % sum;

}



const Velocity& Damp4t10dNotrans::getVelocity() const {
  return *vel;
}

void Damp4t10dNotrans::stepBackward(float* p0, float* p1) const {
  fd4t10s_zjh_2d_vnotrans(p0, p1, &vel->dat[0], vel->nx, vel->nz, dx, dt);
}

void Damp4t10dNotrans::addSource(float* p, const float* source,
    const ShotPosition& pos) const
{
  manipSource(p, source, pos, std::plus<float>());
}

void Damp4t10dNotrans::subSource(float* p, const float* source,
    const ShotPosition& pos) const {
  manipSource(p, source, pos, std::minus<float>());
}

void Damp4t10dNotrans::manipSource(float* p, const float* source,
    const ShotPosition& pos, boost::function2<float, float, float> op) const {
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + nb + EXFDBNDRYLEN;
    int sz = pos.getz(is) + EXFDBNDRYLEN;
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}

void Damp4t10dNotrans::maskGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int bx = nb + EXFDBNDRYLEN;
  int bz = nb + EXFDBNDRYLEN;
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      if (ix < bx || ix >= nxpad - bx || iz >= nzpad - bz) {
        grad[ix  * nzpad + iz] = 0.f;
      }
    }
  }
}

void Damp4t10dNotrans::refillBoundary(float* gradient) const {
  int nz = vel->nz;
  int nx = vel->nx;
  int bx = nb + EXFDBNDRYLEN;
  int bz = nb + EXFDBNDRYLEN;

  const float *srcx0 = gradient + bx * nz;
  const float *srcxn = gradient + (nx - bx - 1) * nz;

  TRACE() << "fill x[0, bx] and x[nx - bx, nx]";
  for (int ix = 0; ix < bx; ix++) {
    float *dstx0 = gradient + ix * nz;
    float *dstxn = gradient + (ix + nx - bx) * nz;
    std::copy(srcx0, srcx0 + nz, dstx0);
    std::copy(srcxn, srcxn + nz, dstxn);
  }

  TRACE() << "fill z[nz - bz, nz]";
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = nz - bx; iz < nz; iz++) {
      int srcIdx = ix * nz + (nz - bz - 1);
      int dstIdx = ix * nz + iz;
      gradient[dstIdx] = gradient[srcIdx];
    }
  }
}

void Damp4t10dNotrans::sfWriteVel(sf_file file) const {
  int nzpad = vel->nz;
  int nxpad = vel->nx;
  int nz = nzpad - 2 * EXFDBNDRYLEN - nb;

  std::vector<float> vv = vel->dat;
  for (int ix = EXFDBNDRYLEN + nb; ix < nxpad - EXFDBNDRYLEN - nb; ix++) {
    sf_floatwrite(&vv[ix * nzpad + EXFDBNDRYLEN], nz, file);
  }
}

void Damp4t10dNotrans::refillVelStencilBndry() {
  Velocity &exvel = getVelocity();
  fillForStencil(exvel, EXFDBNDRYLEN);
}

void Damp4t10dNotrans::removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const {
  int half_len = t_width / dt;
  int sx = allSrcPos.getx(0) + EXFDBNDRYLEN + nb;
  int sz = allSrcPos.getz(0) + EXFDBNDRYLEN;
  int gz = allGeoPos.getz(0) + EXFDBNDRYLEN; // better to assume all receivers are located at the same depth

  float vel_average = 0.0;
  int gmin = (sz < gz) ? sz : gz;
  int gmax = (sz > gz) ? sz : gz;

//  printf("dt %f, half_len %d, sx %d, selav %d, gelav %d\n", dt, half_len, sx, sz, gz);
//  printf("gmin %d, gmax %d\n", gmin, gmax);

  const std::vector<float> &vv = this->vel->dat;
  int nx = this->vel->nx;
  int nz = this->vel->nz;
  for (int i = 0; i < nx; i ++) {
    for (int k = gmin; k <= gmax; k ++) {
      vel_average += vv[i * nz + k];
    }
  }
  vel_average /= nx * (gmax - gmin + 1);
//  printf("vel_average: %.20f\n", vel_average);


  int ng = allGeoPos.ns;
  std::vector<float> trans(nt * ng);
  matrix_transpose(&data[0], &trans[0], ng, nt);

  for (int itr = 0; itr < ng; itr ++) {
    int gx = allGeoPos.getx(itr) + EXFDBNDRYLEN + nb;
    int gz = allGeoPos.getz(itr) + EXFDBNDRYLEN;

    float dist = (gx-sx)*(gx-sx) + (gz-sz)*(gz-sz);
    int t = (int)sqrt(dist * vel_average);
    int start = t;
    int end = ((t + 2 * half_len) > nt) ? nt : (t + 2 * half_len);

    for (int j = start; j < end; j ++) {
      trans[itr * nt + j] = 0.f;
    }
  }

  matrix_transpose(&trans[0], &data[0], nt, ng);
}

Damp4t10dNotrans::Damp4t10dNotrans(const ShotPosition& _allSrcPos, const ShotPosition& _allGeoPos,
    float _dt, float _dx, float _fm, int _nb, int _nt) :
      vel(NULL), allSrcPos(&_allSrcPos), allGeoPos(&_allGeoPos),
      dt(_dt), dx(_dx), fm(_fm), nb(_nb), nt(_nt)
{
  int totalNb = nb + EXFDBNDRYLEN;
  bndr.resize(totalNb);
  initbndr(bndr, bndr.size());
}

void Damp4t10dNotrans::addEncodedSource(float* p, const float* encsrc) const {
  this->addSource(p, encsrc, *this->allSrcPos);
}

void Damp4t10dNotrans::subEncodedSource(float* p, const float* source) const {
  this->subSource(p, source, *this->allSrcPos);
}

void Damp4t10dNotrans::recordSeis(float* seis_it, const float* p) const {
  this->recordSeis(seis_it, p, *this->allGeoPos);
}

void Damp4t10dNotrans::removeDirectArrival(float* data) const {
  float t_width = 1.5 / fm;
  this->removeDirectArrival(*this->allSrcPos, *this->allGeoPos, data, nt, t_width);
}

void Damp4t10dNotrans::addSource(float* p, const float* source, int is) const {

}

int Damp4t10dNotrans::getns() const {
  return allSrcPos->ns;
}

int Damp4t10dNotrans::getng() const {
  return allGeoPos->ns;
}

float Damp4t10dNotrans::getdt() const {
  return dt;
}

float Damp4t10dNotrans::getdx() const {
  return dx;
}

int Damp4t10dNotrans::getnt() const {
  return nt;
}

const ShotPosition& Damp4t10dNotrans::getAllSrcPos() const {
  return *allSrcPos;
}

const ShotPosition& Damp4t10dNotrans::getAllGeoPos() const {
  return *allGeoPos;
}

int Damp4t10dNotrans::getnx() const {
  return vel->nx;
}

int Damp4t10dNotrans::getnz() const {
  return vel->nz;
}

Velocity& Damp4t10dNotrans::getVelocity() {
  return *const_cast<Velocity *>(vel);
}

std::vector<float> Damp4t10dNotrans::initBndryVector(int nt) const {
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

void Damp4t10dNotrans::writeBndry(float* _bndr, const float* p, int it) const {
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

void Damp4t10dNotrans::readBndry(const float* _bndr, float* p, int it) const {
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

