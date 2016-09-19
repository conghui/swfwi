/*
 * damp4t10d.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include <cmath>
#include <functional>
#include "damp4t10d.h"
#include "logger.h"
#include "sum.h"
#include "sfutil.h"
#include "common.h"

extern "C" {
#include <rsf.h>
#include "fd4t10s-damp-zjh.h"
#include "fd4t10s-zjh.h"
}


static void initbndr(std::vector<float> &bndr, int nb) {
  for(int ib=0;ib<nb;ib++){
    float tmp=0.015*(nb-ib);
    bndr[ib]=std::exp(-tmp*tmp);
  }
}

static void applySponge(float* p, const float *bndr, int nx, int nz, int nb) {
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

static void transvel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    vel[i] = (dx * dx) / (vel[i] * vel[i] * dt * dt);
  }
}

static void recoverVel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    vel[i] = std::sqrt(dx*dx / (dt*dt*vel[i]));
  }
}

Velocity Damp4t10d::expandDomain(const Velocity& _vel) {
  // expand for boundary, free surface
  int nb = bx0 - EXFDBNDRYLEN;
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + nb);
  expandBndry(exvelForBndry, _vel, nb);

  transvel(exvelForBndry.dat, dx, dt);

  // expand for stencil
  Velocity ret(exvelForBndry.nx+2*EXFDBNDRYLEN, exvelForBndry.nz+2*EXFDBNDRYLEN);
  expandForStencil(ret, exvelForBndry, EXFDBNDRYLEN);

  return ret;
}

void Damp4t10d::stepForward(float* p0, float* p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);

  fd4t10s_damp_zjh_2d_vtrans(p0, p1, &vel->dat[0], &u2[0], vel->nx, vel->nz, bx0);
}

void Damp4t10d::bindVelocity(const Velocity& _vel) {
  this->vel = &_vel;
}

void Damp4t10d::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) const {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + bx0;
    int gz = geoPos.getz(ig) + bz0;
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }
//  DEBUG() << format("sum %.20f") % sum;

}



const Velocity& Damp4t10d::getVelocity() const {
  return *vel;
}

void Damp4t10d::stepBackward(float* p0, float* p1) const {
  static std::vector<float> u2(vel->nx * vel->nz, 0);
  fd4t10s_zjh_2d_vtrans(p0, p1, &vel->dat[0], &u2[0], vel->nx, vel->nz);
}

void Damp4t10d::addSource(float* p, const float* source,
    const ShotPosition& pos) const
{
  manipSource(p, source, pos, std::plus<float>());
}

void Damp4t10d::subSource(float* p, const float* source,
    const ShotPosition& pos) const {
  manipSource(p, source, pos, std::minus<float>());
}

void Damp4t10d::manipSource(float* p, const float* source,
    const ShotPosition& pos, boost::function2<float, float, float> op) const {
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + bx0;
    int sz = pos.getz(is) + bz0;
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}

void Damp4t10d::maskGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      if (ix < bx0 || iz < bz0 || ix >= nxpad - bxn || iz >= nzpad - bzn) {
        grad[ix  * nzpad + iz] = 0.f;
      }
    }
  }
}

void Damp4t10d::refillBoundary(float* gradient) const {
  int nzpad = vel->nz;
  int nxpad = vel->nx;

  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < bz0; iz++) {
      gradient[ix * nzpad + iz] = gradient[ix * nzpad + bz0];           // top
    }
    for (int iz = nzpad - bzn; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[ix * nzpad + nzpad - bzn - 1];  // bottom
    }
  }

  for (int ix = 0; ix < bx0; ix++) {                              // left
    for (int iz = 0; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[bx0 * nzpad + iz];
    }
  }

  for (int ix = nxpad - bxn; ix < nxpad; ix++) {                        // right
    for (int iz = 0; iz < nzpad; iz++) {
      gradient[ix * nzpad + iz] = gradient[(nxpad - bxn - 1) * nzpad + iz];
    }
  }
}

void Damp4t10d::sfWriteVel(const std::vector<float> &exvel, sf_file file) const {
  assert(exvel.size() == vel->dat.size());
  int nzpad = vel->nz;
  int nxpad = vel->nx;
  int nz = nzpad - bz0 - bzn;

  std::vector<float> vv = exvel;
  recoverVel(vv, dx, dt);
  for (int ix = bx0; ix < nxpad - bxn; ix++) {
    sf_floatwrite(&vv[ix * nzpad + EXFDBNDRYLEN], nz, file);
  }
}

void Damp4t10d::refillVelStencilBndry() {
  Velocity &exvel = getVelocity();
  fillForStencil(exvel, EXFDBNDRYLEN);
}

void Damp4t10d::FwiForwardModeling(const std::vector<float>& encSrc,
    std::vector<float>& dcal, int shot_id) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);
  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);

  for(int it=0; it<nt; it++) {
    addSource(&p1[0], &encSrc[it], curSrcPos);
    stepForward(&p0[0], &p1[0]);
    std::swap(p1, p0);
    recordSeis(&dcal[it*ng], &p0[0]);
  }
}

void Damp4t10d::EssForwardModeling(const std::vector<float>& encSrc,
    std::vector<float>& dcal) const {
  int nx = getnx();
  int nz = getnz();
  int ns = getns();
  int ng = getng();

  std::vector<float> p0(nz * nx, 0);
  std::vector<float> p1(nz * nx, 0);

  for(int it=0; it<nt; it++) {
    addEncodedSource(&p1[0], &encSrc[it * ns]);
    stepForward(&p0[0], &p1[0]);
    std::swap(p1, p0);
    recordSeis(&dcal[it*ng], &p0[0]);
  }
}

void Damp4t10d::scaleGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int nz = nzpad - bz0 - bzn;

  for (int ix = bx0; ix < nxpad - bxn; ix++) {
    for (int iz = 1; iz < nz; iz++) {
      grad[ix*nzpad + iz+bz0] *= std::sqrt(static_cast<float>(iz));
    }
  }
}

void Damp4t10d::removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const {
  int half_len = t_width / dt;
  int sx = allSrcPos.getx(0) + bx0;
  int sz = allSrcPos.getz(0) + bz0;
  int gz = allGeoPos.getz(0) + bz0; // better to assume all receivers are located at the same depth

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
    int gx = allGeoPos.getx(itr) + bx0;
    int gz = allGeoPos.getz(itr) + bz0;

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

Damp4t10d::Damp4t10d(const ShotPosition& _allSrcPos, const ShotPosition& _allGeoPos,
    float _dt, float _dx, float _fm, int _nb, int _nt) :
      vel(NULL), allSrcPos(&_allSrcPos), allGeoPos(&_allGeoPos),
      dt(_dt), dx(_dx), fm(_fm),  nt(_nt)
{
  bz0 = EXFDBNDRYLEN;
  bx0 = bxn = bzn = _nb + EXFDBNDRYLEN;
  bndr.resize(bx0);
  initbndr(bndr, bndr.size());
}

void Damp4t10d::addEncodedSource(float* p, const float* encsrc) const {
  this->addSource(p, encsrc, *this->allSrcPos);
}

void Damp4t10d::subEncodedSource(float* p, const float* source) const {
  this->subSource(p, source, *this->allSrcPos);
}

void Damp4t10d::recordSeis(float* seis_it, const float* p) const {
  this->recordSeis(seis_it, p, *this->allGeoPos);
}

void Damp4t10d::fwiRemoveDirectArrival(float* data, int shot_id) const {
  float t_width = 1.5 / fm;
  ShotPosition curSrcPos = allSrcPos->clipRange(shot_id, shot_id);
  this->removeDirectArrival(curSrcPos, *this->allGeoPos, data, nt, t_width);
}

void Damp4t10d::removeDirectArrival(float* data) const {
  float t_width = 1.5 / fm;
  this->removeDirectArrival(*this->allSrcPos, *this->allGeoPos, data, nt, t_width);
}

void Damp4t10d::addSource(float* p, const float* source, int is) const {

}

int Damp4t10d::getns() const {
  return allSrcPos->ns;
}

int Damp4t10d::getng() const {
  return allGeoPos->ns;
}

float Damp4t10d::getdt() const {
  return dt;
}

float Damp4t10d::getdx() const {
  return dx;
}

int Damp4t10d::getnt() const {
  return nt;
}

const ShotPosition& Damp4t10d::getAllSrcPos() const {
  return *allSrcPos;
}

const ShotPosition& Damp4t10d::getAllGeoPos() const {
  return *allGeoPos;
}

int Damp4t10d::getnx() const {
  assert(vel != NULL);
  return vel->nx;
}

int Damp4t10d::getnz() const {
  assert(vel != NULL);
  return vel->nz;
}

Velocity& Damp4t10d::getVelocity() {
  return *const_cast<Velocity *>(vel);
}

std::vector<float> Damp4t10d::initBndryVector(int nt) const {
  if (vel == NULL) {
    ERROR() << __PRETTY_FUNCTION__ << ": you should bind velocity first";
    exit(1);
  }
  int nxpad = vel->nx;
  int nzpad = vel->nz;
//  int nx = nxpad - 2*nb;
//  int nz = nzpad - nb;

//  int nx = nxpad - bx0 - bxn;
//  int nz = nzpad - bz0 - bzn;
//
//  bndrLen = FDLEN + 1;
//  bndrSize = bndrLen * (
//             nx +   /* bottom */
//             2 * nz /* left + right */
//             );

  bndrWidth = 6;
  int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
  int nz = nzpad - bz0 - bzn;

  bndrSize = bndrWidth * (
      nx +   /* bottom */
      2 * nz /* left + right */
  );

  return std::vector<float>(nt*bndrSize, 0);
}

void Damp4t10d::writeBndry(float* _bndr, const float* p, int it) const {
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
     *    *******************
     *    *******************
     *
     */
    int nxpad = vel->nx;
    int nzpad = vel->nz;

    int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
    int nz = nzpad - bz0 - bzn;

    float *bndr = &_bndr[it * bndrSize];

    for (int ix = 0; ix < nx; ix++) {
      for(int iz = 0; iz < bndrWidth; iz++) {
        bndr[iz + bndrWidth*ix] = p[(ix+bx0-bndrWidth)*nzpad + (nzpad - bzn + iz)]; // bottom
      }
    }

    for (int iz = 0; iz < nz; iz++) {
      for(int ix=0; ix < bndrWidth; ix++) {
        bndr[bndrWidth*nx+iz+nz*ix]         = p[(bx0-bndrWidth + ix)*nzpad + (bz0 + iz)];   // left
        bndr[bndrWidth*nx+iz+nz*(ix+bndrWidth)] = p[(nxpad - bxn + ix)*nzpad + (bz0 + iz)];  // right
      }
    }
}

void Damp4t10d::readBndry(const float* _bndr, float* p, int it) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;

  int nx = nxpad - (bx0 - bndrWidth + bxn - bndrWidth);
  int nz = nzpad - bz0 - bzn;
  const float *bndr = &_bndr[it * bndrSize];

  for (int ix = 0; ix < nx; ix++) {
    for(int iz = 0; iz < bndrWidth; iz++) {
      p[(ix+bx0-bndrWidth)*nzpad + (nzpad - bzn + iz)] = bndr[iz + bndrWidth*ix]; // bottom
    }
  }

  for (int iz = 0; iz < nz; iz++) {
    for(int ix=0; ix < bndrWidth; ix++) {
      p[(bx0-bndrWidth + ix)*nzpad + (bz0 + iz)] = bndr[bndrWidth*nx+iz+nz*ix];   // left
      p[(nxpad - bxn + ix)*nzpad + (bz0 + iz)] = bndr[bndrWidth*nx+iz+nz*(ix+bndrWidth)];  // right
    }
  }
}

