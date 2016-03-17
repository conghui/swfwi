/*
 * ReguFactor.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: Conghui He (heconghui@gmail.com)
 */

#include "ReguFactor.h"

#include <algorithm>
#include <functional>

ReguFactor::ReguFactor(const float *vel, int nx, int nz, float lambdaX, float lambdaZ) :
  mModel(vel), mNx(nx), mNz(nz), mLambdaX(lambdaX), mLambdaZ(lambdaZ),
  mReguTermValue(0), mReguGradVector(nx * nz) {
}

float ReguFactor::getReguTerm() {
  float Wx2 = getWx2();
  float Wz2 = getWz2();

  /// update reguTermValue
  mReguTermValue = mLambdaX * Wx2 + mLambdaZ * Wz2;

  return mReguTermValue;
}

float ReguFactor::getWx2() const {
  /// calculate Wx2
  float Wx2 = 0.0f;

  const int nx = mNx;
  const int nz = mNz;
  const float *p = mModel;

  for (int ix = 0; ix < nx - 1; ix++) { /// nx - 1
    for (int iz = 0; iz < nz; iz++) {
      int idx     = (ix + 0) * nz + iz;
      int idx_xp1 = (ix + 1) * nz + iz;

      float r = p[idx_xp1] - p[idx];
      Wx2    += r * r;
    }
  }

  return Wx2;
}

float ReguFactor::getWz2() const {
  /// calculate Wz2
  float Wz2 = 0.0f;

  const int nx = mNx;
  const int nz = mNz;
  const float *p = mModel;

  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz - 1; iz++) { /// nz - 1
      int idx     = ix * nz + iz + 0;
      int idx_zp1 = ix * nz + iz + 1;

      float r = p[idx_zp1] - p[idx];
      Wz2    += r * r;
    }
  }

  return Wz2;
}

const float *ReguFactor::getReguGradient() {
  const int nx = mNx;
  const int nz = mNz;
  const float *p = mModel;

  std::vector<float> xTerm(nx * nz, 0);
  std::vector<float> zTerm(nx * nz, 0);

  for (int ix = 0; ix < nx - 1; ix++) { /// nx - 1
    for (int iz = 0; iz < nz; iz++) {
      int idx     = (ix + 0) * nz + iz;
      int idx_xp1 = (ix + 1) * nz + iz;
      int idx_xm1 = (ix - 1) * nz + iz;

      if (ix == 0) {
        xTerm[idx] = 2 * mLambdaX * (p[idx] - p[idx_xp1]);
      } else {
        xTerm[idx] = 2 * mLambdaX * (2 * p[idx] - p[idx_xp1] - p[idx_xm1]);
      }
    }
  }

  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz - 1; iz++) { /// nz - 1
      int idx     = ix * nz + iz + 0;
      int idx_zp1 = ix * nz + iz + 1;
      int idx_zm1 = ix * nz + iz - 1;

      if (iz == 0) {
        xTerm[idx] = 2 * mLambdaZ * (p[idx] - p[idx_zp1]);
      } else {
        xTerm[idx] = 2 * mLambdaZ * (2 * p[idx] - p[idx_zp1] - p[idx_zm1]);
      }
    }
  }

  std::transform(xTerm.begin(), xTerm.end(), zTerm.begin(), mReguGradVector.begin(), std::plus<float>());

  return &mReguGradVector[0];
}
