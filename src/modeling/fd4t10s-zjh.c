/*
 * fd4t10s-zjh.cpp
 *
 *  Created on: Mar 5, 2016
 *      Author: rice
 */

#include "fd4t10s-zjh.h"

void fd4t10s_zjh_2d_vtrans(float *prev_wave, const float *curr_wave, const float *vel, float *u2, int nx, int nz) {
  float a[6];

  const int d = 5;
  int ix, iz;

  /// Zhang, Jinhai's method
  a[0] = +1.53400796;
  a[1] = +1.78858721;
  a[2] = -0.31660756;
  a[3] = +0.07612173;
  a[4] = -0.01626042;
  a[5] = +0.00216736;

  #pragma omp parallel for default(shared) private(ix, iz)
  for (ix = d; ix < nx - d; ix ++) {
    for (iz = d; iz < nz - d; iz ++) {
      int curPos = ix * nz + iz;
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

  #pragma omp parallel for default(shared) private(ix, iz)
  for (ix = d + 1; ix < nx - d - 1; ix++) { /// the range of ix is different from that in previous for loop
    for (iz = d + 1; iz < nz - d - 1; iz++) { /// be careful of the range of iz
      int curPos = ix * nz + iz;
      float curvel = vel[curPos];

      prev_wave[curPos] = 2. * curr_wave[curPos] - prev_wave[curPos]  +
                          (1.0f / curvel) * u2[curPos] + /// 2nd order
                          1.0f / 12 * (1.0f / curvel) * (1.0f / curvel) *
                          (u2[curPos - 1] + u2[curPos + 1] + u2[curPos - nz] + u2[curPos + nz] - 4 * u2[curPos]); /// 4th order
    }
  }

}
