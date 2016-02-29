/*
 * damp4t10d.cpp
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#include "damp4t10d.h"

Damp4t10d::Damp4t10d() {
  // TODO Auto-generated constructor stub

}


//static void fd4t10s_bd1_2d(dim3d_t dim, const float *vel, float *prev_wave, float *curr_wave) {
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
