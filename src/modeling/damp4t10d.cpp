/*
 * damp4t10d.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include "damp4t10d.h"
#include "logger.h"
#include "sum.h"
#include "fd4t10s-damp-zjh.h"

extern "C" {
#include <rsf.h>
}

static void expandForStencil(Velocity &exvel, const Velocity &v0, int halo) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nxpad = nx + 2 * halo;
  int nzpad = nz + 2 * halo;

  std::vector<float> &vel_e = exvel.dat;
  const std::vector<float> &vel = v0.dat;

  //copy the vel into vel_e
  for (int ix = halo; ix < nx + halo; ix++) {
    for (int iz = halo; iz < nz + halo; iz++) {
      vel_e[ix * nzpad + iz] = vel[(ix - halo) * nz +  (iz - halo)];
    }
  }

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
      b[(nb + ix) * nzpad + iz] = a[ix * nx + iz];
    }
  }

  /// boundary
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + iz] = b[ix * nzpad + nb];                              /* top */
      b[ix * nzpad + (nzpad - iz - 1)] = b[ix * nzpad + (nzpad - nb - 1)];  /* bottom*/
    }
  }

  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      b[ix * nzpad + iz] = b[nb * nzpad + iz];                              /* left */
      b[(nxpad - ix - 1) * nzpad + iz] = b[(nxpad - nb - 1) * nzpad + iz];  /* right */
    }
  }
}

Damp4t10d::Damp4t10d(float _dt, float _dx, int _nb) :
  vel(NULL), dt(_dt), dx(_dx), nb(_nb)
{

}

Velocity Damp4t10d::expandDomain(const Velocity& _vel) {
  Velocity vv = _vel;

  // expand for boundary
  Velocity exvelForBndry(vv.nx + 2 * nb, vv.nz + nb);
  expandBndry(exvelForBndry, vv, nb);

  // expand for stencil
  Velocity exvelForStencil(exvelForBndry.nx+2*FDLEN, exvelForBndry.nz+2*FDLEN);
  expandForStencil(exvelForStencil, vv, FDLEN);

  DEBUG() << format("sum exvel %.20f") % sum(exvelForStencil.dat);

  return exvelForStencil;
}

void Damp4t10d::stepForward(float* p0, float* p1) {
  fd4t10s_damp_zjh_2d(p0, p1, &vel->dat[0], vel->nx, vel->nz, nb, dx, dt);
}

void Damp4t10d::bindVelocity(const Velocity& _vel) {
  this->vel = &_vel;
}

void Damp4t10d::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) {

  int ng = geoPos.ns;
  int nzpad = vel->nz;

//  DEBUG() << format("ng %d") % ng;
//  float sum = 0;
  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig) + nb + FDLEN;
    int gz = geoPos.getz(ig) + FDLEN;
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
//    DEBUG() << format("ig %d, idx %d, v %.20f") % ig % idx % seis_it[ig];
  }
//  DEBUG() << format("sum %.20f") % sum;

}

void Damp4t10d::addSource(float* p, const float* source,
    const ShotPosition& pos)
{
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is) + nb + FDLEN;
    int sz = pos.getz(is) + FDLEN;
    p[sx * nzpad + sz] += source[is];
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}
