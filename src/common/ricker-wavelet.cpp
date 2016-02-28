/*
 * ricker-wavelet.cpp
 *
 *  Created on: Feb 24, 2016
 *      Author: rice
 */

extern "C" {
#include <rsf.h>
}
#include <cmath>

#include "ricker-wavelet.h"

void rickerWavelet(float *wlt, int nt, float fm, float dt) {
  for (int it = 0; it < nt; it++) {
    float tmp = M_PI * fm * (it * dt - 1.0 / fm);
    tmp *= tmp;
    wlt[it] = (1.0 - 2.0 * tmp) * std::exp(-tmp);
  }
}
