/*
 * damp4t10d.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include <functional>
#include "damp4t10d.h"
#include "logger.h"
#include "sum.h"
#include "fd4t10s-damp-zjh.h"
#include "fd4t10s-zjh.h"
#include "sfutil.h"
#include "common.h"

extern "C" {
#include <rsf.h>
}

static void expandForStencil(Velocity &exvel, const Velocity &v0, int halo) {
  int nx = v0.nx;
  int nz = v0.nz;
  int nxpad = nx + 2 * halo;
  int nzpad = nz + 2 * halo;

  const std::vector<float> &vel = v0.dat;
  std::vector<float> &vel_e = exvel.dat;

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


Damp4t10d::Damp4t10d(float _dt, float _dx, int _nb) :
  vel(NULL), dt(_dt), dx(_dx), nb(_nb)
{

}

static void transvel(std::vector<float> &vel, float dx, float dt) {
  for (size_t i = 0; i < vel.size(); i ++) {
    vel[i] = (dx * dx) / (vel[i] * vel[i] * dt * dt);
  }
}

Velocity Damp4t10d::expandDomain(const Velocity& _vel) {
  // expand for boundary, free surface
  Velocity exvelForBndry(_vel.nx + 2 * nb, _vel.nz + nb);
  expandBndry(exvelForBndry, _vel, nb);
  sfFloatWrite2d("000vel.rsf", &exvelForBndry.dat[0], exvelForBndry.nz, exvelForBndry.nx);

  transvel(exvelForBndry.dat, dx, dt);

  // expand for stencil
  Velocity ret(exvelForBndry.nx+2*FDLEN, exvelForBndry.nz+2*FDLEN);
  expandForStencil(ret, exvelForBndry, FDLEN);

  return ret;
}

void Damp4t10d::stepForward(float* p0, float* p1) const {
  fd4t10s_damp_zjh_2d(p0, p1, &vel->dat[0], vel->nx, vel->nz, nb + FDLEN, dx, dt);
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
    int gx = geoPos.getx(ig) + nb + FDLEN;
    int gz = geoPos.getz(ig) + FDLEN;
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
  fd4t10s_zjh_2d(p0, p1, &vel->dat[0], vel->nx, vel->nz, nb + FDLEN, dx, dt);
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
    int sx = pos.getx(is) + nb + FDLEN;
    int sz = pos.getz(is) + FDLEN;
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
//    DEBUG() << format("sx %d, sz %d, source[%d] %f") % sx % sz % is % source[is];
  }
}

void Damp4t10d::maskGradient(float* grad) const {
  int nxpad = vel->nx;
  int nzpad = vel->nz;
  int bx = nb + FDLEN;
  int bz = nb + FDLEN;
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      if (ix < bx || ix >= nxpad - bx || iz >= nzpad - bz) {
        grad[ix  * nzpad + iz] = 0.f;
      }
    }
  }
}

void Damp4t10d::refillBoundary(float* gradient) const {
  int nz = vel->nz;
  int nx = vel->nx;
  int bx = nb + FDLEN;
  int bz = nb + FDLEN;

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

void Damp4t10d::sfWriteVel(sf_file file) const {
  int nzpad = vel->nz;
  int nxpad = vel->nx;
  int nz = nzpad - 2 * FDLEN - nb;

  for (int ix = FDLEN + nb; ix < nxpad - FDLEN - nb; ix++) {
    sf_floatwrite(const_cast<float *>(&vel->dat[ix * nzpad + FDLEN]), nz, file);
  }
}

void Damp4t10d::removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const {
  int half_len = t_width / dt;
  int sx = allSrcPos.getx(0) + FDLEN + nb;
  int sz = allSrcPos.getz(0) + FDLEN;
  int gz = allGeoPos.getz(0) + FDLEN; // better to assume all receivers are located at the same depth

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
    int gx = allGeoPos.getx(itr) + FDLEN + nb;
    int gz = allGeoPos.getz(itr) + FDLEN;

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
