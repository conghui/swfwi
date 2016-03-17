/*
 * enquistabc2d.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: rice
 */

#include <cstdlib>
#include <functional>
#include "enquistabc2d.h"

EnquistAbc2d::EnquistAbc2d(float dt, float dx, float dz) :
  IModeling(), dt(dt), dx(dx), dz(dz), dtx(dt/dx), dtz(dt/dz)
{

}

void EnquistAbc2d::stepForward(const float *p0, const float *p1, float *p2) const {
  int nx = vel->nx;
  int nz = vel->nz;
  const std::vector<float> &vv = vel->dat;

#pragma omp parallel for
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      float v1 = vv[ix * nz + iz] * dtz;
      v1 = v1 * v1;
      float v2 = vv[ix * nz + iz] * dtx;
      v2 = v2 * v2;
      float diff1 = -2.0 * p1[ix * nz + iz];
      float diff2 = diff1;

      diff1 += (iz - 1 >= 0) ? p1[ix * nz + (iz - 1)] : 0.0;
      diff1 += (iz + 1 < nz) ? p1[ix * nz + (iz + 1)] : 0.0;
      diff2 += (ix - 1 >= 0) ? p1[(ix - 1) * nz + iz] : 0.0;
      diff2 += (ix + 1 < nx) ? p1[(ix + 1) * nz + iz] : 0.0;
      diff1 *= v1;
      diff2 *= v2;
      p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
    }
  }

#pragma omp parallel for
  for (int ix = 1; ix < nx - 1; ix++) {
    /* top boundary */
    /*
    iz=0;
    diff1= (p1[ix][iz+1]-p1[ix][iz])-
    (p0[ix][iz+1]-p0[ix][iz]);
    diff2= c21*(p1[ix-1][iz]+p1[ix+1][iz]) +
    c22*(p1[ix-2][iz]+p1[ix+2][iz]) +
    c20*p1[ix][iz];
    diff1*=sqrtf(vv[ix][iz])/dz;
    diff2*=vv[ix][iz]/(2.0*dx*dx);
    p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
     */
    /* bottom boundary */
    int iz = nz - 1;
    float v1 = vv[ix * nz + iz] * dtz;
    float v2 = vv[ix * nz + iz] * dtx;
    float diff1 = -(p1[ix * nz + iz] - p1[ix * nz + (iz - 1)]) + (p0[ix * nz + iz] - p0[ix * nz + iz - 1]);
    float diff2 = p1[(ix - 1) * nz + iz] - 2.0 * p1[ix * nz + iz] + p1[(ix + 1) * nz + iz];
    diff1 *= v1;
    diff2 *= 0.5 * v2 * v2;
    p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
  }

#pragma omp parallel for
  for (int iz = 1; iz < nz - 1; iz++) {
    /* left boundary */
    int ix = 0;
    float v1 = vv[ix * nz + iz] * dtz;
    float v2 = vv[ix * nz + iz] * dtx;
    float diff1 = p1[ix * nz + (iz - 1)] - 2.0 * p1[ix * nz + iz] + p1[ix * nz + (iz + 1)];
    float diff2 = (p1[(ix + 1) * nz + iz] - p1[ix * nz + iz]) - (p0[(ix + 1) * nz + iz] - p0[ix * nz + iz]);
    diff1 *= 0.5 * v1 * v1;
    diff2 *= v2;
    p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;

    /* right boundary */
    ix = nx - 1;
    v1 = vv[ix * nz + iz] * dtz;
    v2 = vv[ix * nz + iz] * dtx;
    diff1 = p1[ix * nz + (iz - 1)] - 2.0 * p1[ix * nz + iz] + p1[ix * nz + (iz + 1)];
    diff2 = -(p1[ix * nz + iz] - p1[(ix - 1) * nz + iz]) + (p0[ix * nz + iz] - p0[(ix - 1) * nz + iz]);
    diff1 *= 0.5 * v1 * v1;
    diff2 *= v2;
    p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
  }
}

void EnquistAbc2d::addSource(float* p, const float* source,
    const ShotPosition& pos) const {

  manipSource(p, source, pos, std::plus<float>());
}

void EnquistAbc2d::recordSeis(float* seis_it, const float* p,
    const ShotPosition& geoPos) const {
  int ng = geoPos.ns;
  int nzpad = vel->nz;

  for (int ig = 0; ig < ng; ig++) {
    int gx = geoPos.getx(ig);
    int gz = geoPos.getz(ig);
    int idx = gx * nzpad + gz;
    seis_it[ig] = p[idx];
  }
}

void EnquistAbc2d::writeBndry(float* bndr, const float* p, int it) const {
  int nx = vel->nx;
  int nz = vel->nz;

  float *bndrNow = &bndr[it * (2 * nz + nx)];
  for (int i = 0; i < nz; i++) {
    bndrNow[i] = p[i];
    bndrNow[i + nz] = p[(nx - 1) * nz + i];
  }
  for (int i = 0; i < nx; i++) {
    bndrNow[i + 2 * nz] = p[i * nz + (nz - 1)];
  }
}

void EnquistAbc2d::readBndry(const float* bndr, float* p, int it) const {
  int nx = vel->nx;
  int nz = vel->nz;

  const float *bndrNow = &bndr[it * (2 * nz + nx)];
  for (int i = 0; i < nz; i++) {
    p[i] = bndrNow[i];
    p[(nx - 1) * nz + i] = bndrNow[i + nz];
  }
  for (int i = 0; i < nx; i++) {
    p[i * nz + (nz - 1)] = bndrNow[i + 2 * nz];
  }
}

void EnquistAbc2d::subSource(float* p, const float* source,
    const ShotPosition& pos) const {

  manipSource(p, source, pos, std::minus<float>());
}

std::vector<float> EnquistAbc2d::initBndryVector(int nt) const {
  int nx = vel->nx;
  int nz = vel->nz;
  return std::vector<float>(nt * (2 * nz + nx), 0); /* boundaries for wavefield reconstruction */
}

void EnquistAbc2d::stepBackward(float* illum, float* lap, const float* p0,
    const float* p1, float* p2) const {
  int nx = vel->nx;
  int nz = vel->nz;
  const std::vector<float> &vv = vel->dat;

#pragma omp parallel for
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      float v1 = vv[ix * nz + iz] * dtz;
      v1 = v1 * v1;
      float v2 = vv[ix * nz + iz] * dtx;
      v2 = v2 * v2;
      float diff1 = -2.0 * p1[ix * nz + iz];
      float diff2 = diff1;
      diff1 += (iz - 1 >= 0) ? p1[ix * nz + (iz - 1)] : 0.0;
      diff1 += (iz + 1 < nz) ? p1[ix * nz + (iz + 1)] : 0.0;
      diff2 += (ix - 1 >= 0) ? p1[(ix - 1) * nz + iz] : 0.0;
      diff2 += (ix + 1 < nx) ? p1[(ix + 1) * nz + iz] : 0.0;
      lap[ix * nz + iz] = diff1 + diff2;
      diff1 *= v1;
      diff2 *= v2;
      p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
      illum[ix * nz + iz] += p1[ix * nz + iz] * p1[ix * nz + iz];
    }
  }
}

Velocity EnquistAbc2d::expandDomain(const Velocity& v0) const {
  return v0;
}

void EnquistAbc2d::manipSource(float* p, const float* source,
    const ShotPosition& pos, boost::function2<float, float, float> op) const {
  int nzpad = vel->nz;

  for (int is = 0; is < pos.ns; is++) {
    int sx = pos.getx(is);
    int sz = pos.getz(is);
//    p[sx * nzpad + sz] += source[is];
    p[sx * nzpad + sz] = op(p[sx * nzpad + sz], source[is]);
  }
}
