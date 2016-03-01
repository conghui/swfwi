/*
 * damp4t10d.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include "damp4t10d.h"
#include "logger.h"
#include "sum.h"

extern "C" {
#include <rsf.h>
}

//static void fd4t10s_bd1_2d(float *prev_wave, const float *curr_wave) {
//  float a[6];
//
//  int d = 5;
//
//  int ix, iz, curPos;
//
//  /// Zhang, Jinhai's method
//  a[0] = +1.53400796;
//  a[1] = +1.78858721;
//  a[2] = -0.31660756;
//  a[3] = +0.07612173;
//  a[4] = -0.01626042;
//  a[5] = +0.00216736;
//
//  float max_delta = 0.05;
//
//  int nx =
//  std::vector<float> u2(dim.nx * dim.nz, 0);
//
//  #pragma omp parallel for schedule(dynamic) private(ix, iz, curPos)
//  for (ix = d; ix < dim.nx - d; ix ++) {
//    for (iz = d; iz < dim.nz - d; iz ++) {
//      curPos = idx_in_cube(dim, ix, 0, iz);
//      u2[curPos] = -4.0 * a[0] * curr_wave[curPos] +
//                   a[1] * (curr_wave[curPos - 1]  +  curr_wave[curPos + 1]  +
//                           curr_wave[curPos - dim.nz]  +  curr_wave[curPos + dim.nz])  +
//                   a[2] * (curr_wave[curPos - 2]  +  curr_wave[curPos + 2]  +
//                           curr_wave[curPos - 2 * dim.nz]  +  curr_wave[curPos + 2 * dim.nz])  +
//                   a[3] * (curr_wave[curPos - 3]  +  curr_wave[curPos + 3]  +
//                           curr_wave[curPos - 3 * dim.nz]  +  curr_wave[curPos + 3 * dim.nz])  +
//                   a[4] * (curr_wave[curPos - 4]  +  curr_wave[curPos + 4]  +
//                           curr_wave[curPos - 4 * dim.nz]  +  curr_wave[curPos + 4 * dim.nz])  +
//                   a[5] * (curr_wave[curPos - 5]  +  curr_wave[curPos + 5]  +
//                           curr_wave[curPos - 5 * dim.nz]  +  curr_wave[curPos + 5 * dim.nz]);
//    }
//  }
//
//  #pragma omp parallel for schedule(dynamic) private(ix, iz, curPos)
//  for (ix = d + 1; ix < dim.nx - d - 1; ix ++) { /// the range of ix is different from that in previous for loop
//    for (iz = d + 1; iz < dim.nz - d - 1; iz ++) { /// be careful of the range of iz
//      float delta;
//      float dist = 0;
//      if (ix >= dim.bx && ix < dim.nx - dim.bx &&
//          iz < dim.nz - dim.bz) {
//        dist = 0;
//      }
//      if (ix < dim.bx) {
//        dist = (float)(dim.bx - ix) / dim.bx;
//      }
//      if (ix >= dim.nx - dim.bx) {
//        dist = (float)(ix - (dim.nx - dim.bx)  +  1) / dim.bx;
//      }
//      if (iz >= dim.nz - dim.bz) {
//        dist = (float)(iz - (dim.nz - dim.bz) + 1) / dim.bz;
//      }
//
//      delta = max_delta * dist * dist;
//
//      curPos = idx_in_cube(dim, ix, 0, iz);
//      prev_wave[curPos] = (2. - 2 * delta + delta * delta) * curr_wave[curPos] - (1 - 2 * delta) * prev_wave[curPos]  +
//                          (1.0f / vel[curPos]) * u2[curPos] + /// 2nd order
//                          1.0f / 12 * (1.0f / vel[curPos]) * (1.0f / vel[curPos]) *
//                          (u2[curPos - 1] + u2[curPos + 1] + u2[curPos - dim.nz] + u2[curPos + dim.nz] - 4 * u2[curPos]); /// 4th order
//    }
//  }
//}


static void numericTrans(Velocity &vel, float dx, float dt) {
  std::vector<float> &vv = vel.dat;
  for (size_t i = 0; i < vv.size(); i ++) {
    vv[i] = (dx * dx) / (vv[i] * vv[i] * dt * dt);
  }
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
  dt(_dt), dx(_dx), nb(_nb)
{

}

Velocity Damp4t10d::transformVelocityForModeling(const Velocity& _vel) {
  Velocity vv = _vel;

  DEBUG() << format("nx %d, nz %d") % vv.nx % vv.nz;
  DEBUG() << format("before transform sum vv %.20f") % sum(vv.dat);
  numericTrans(vv, dx, dt);

  DEBUG() << format("sum vv %.20f") % sum(vv.dat);

  // expand for boundary
  Velocity exvelForBndry(vv.nx + 2 * nb, vv.nz + nb);
  expandBndry(exvelForBndry, vv, nb);

  Velocity exvelForStencil(exvelForBndry.nx+2*FDLEN, exvelForBndry.nz+2*FDLEN);
  expandForStencil(exvelForStencil, vv, FDLEN);

  sf_file f = sf_output("exvel.rsf");
  sf_putint(f, "n1", exvelForStencil.nz);
  sf_putint(f, "n2", exvelForStencil.nx);
  sf_putfloat(f, "d1", dx);
  sf_putfloat(f, "d2", dx);
  sf_putfloat(f, "o1", 0);
  sf_putfloat(f, "o2", 0);
  sf_floatwrite(&exvelForStencil.dat[0], exvelForStencil.nx * exvelForStencil.nz, f);

  DEBUG() << format("sum exvel %.20f") % sum(exvelForStencil.dat);

  return vv;
}

void Damp4t10d::stepForward(float* p0, float* p1) {
  fd4t10s_bd1_2d(p0, p1);
}

void Damp4t10d::setVelocity(const Velocity& _vel) {
  this->vel = &_vel;
}

void Damp4t10d::fd4t10s_bd1_2d(float* prev_wave, const float* curr_wave) {
  float a[6];

  int d = 5;

  int ix, iz, curPos;

  /// Zhang, Jinhai's method
  a[0] = +1.53400796;
  a[1] = +1.78858721;
  a[2] = -0.31660756;
  a[3] = +0.07612173;
  a[4] = -0.01626042;
  a[5] = +0.00216736;

  float max_delta = 0.05;

  int nx = vel->nx;
  int nz = vel->nz;
  std::vector<float> u2(nx * nz, 0);

  #pragma omp parallel for schedule(dynamic) private(ix, iz, curPos)
  for (ix = d; ix < nx - d; ix ++) {
    for (iz = d; iz < nz - d; iz ++) {
//      curPos = idx_in_cube(dim, ix, 0, iz);
      curPos = ix * nz + iz;
      u2[curPos] = -4.0 * a[0] * curr_wave[curPos] +
                   a[1] * (curr_wave[curPos - 1]  +  curr_wave[curPos + 1]  +
                           curr_wave[curPos - nz]  +  curr_wave[curPos + nz])  +
                   a[2] * (curr_wave[curPos - 2]  +  curr_wave[curPos + 2]  +
                           curr_wave[curPos - 2 * nz]  +  curr_wave[curPos + 2 * nz])  +
                   a[3] * (curr_wave[curPos - 3]  +  curr_wave[curPos + 3]  +
                           curr_wave[curPos - 3 * nz]  +  curr_wave[curPos + 3 * nz])  +
                   a[4] * (curr_wave[curPos - 4]  +  curr_wave[curPos + 4]  +
                           curr_wave[curPos - 4 * nz]  +  curr_wave[curPos + 4 * nz])  +
                   a[5] * (curr_wave[curPos - 5]  +  curr_wave[curPos + 5]  +
                           curr_wave[curPos - 5 * nz]  +  curr_wave[curPos + 5 * nz]);
    }
  }

  #pragma omp parallel for schedule(dynamic) private(ix, iz, curPos)
  for (ix = d + 1; ix < nx - d - 1; ix ++) { /// the range of ix is different from that in previous for loop
    for (iz = d + 1; iz < nz - d - 1; iz ++) { /// be careful of the range of iz
      float delta;
      float dist = 0;
      if (ix >= nb && ix < nx - nb && iz < nz - nb) {
        dist = 0;
      }
      if (ix < nb) {
        dist = (float)(nb - ix) / nb;
      }
      if (ix >= nx - nb) {
        dist = (float)(ix - (nx - nb)  +  1) / nb;
      }
      if (iz >= nz - nb) {
        dist = (float)(iz - (nz - nb) + 1) / nb;
      }

      delta = max_delta * dist * dist;

//      curPos = idx_in_cube(dim, ix, 0, iz);
      curPos = ix * nz + iz;
      float curvel = vel->dat[curPos];
      prev_wave[curPos] = (2. - 2 * delta + delta * delta) * curr_wave[curPos] - (1 - 2 * delta) * prev_wave[curPos]  +
                          (1.0f / curvel) * u2[curPos] + /// 2nd order
                          1.0f / 12 * (1.0f / curvel) * (1.0f / curvel) *
                          (u2[curPos - 1] + u2[curPos + 1] + u2[curPos - nz] + u2[curPos + nz] - 4 * u2[curPos]); /// 4th order
    }
  }
}
