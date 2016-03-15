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

void rickerWavelet(float *wlt, int nt, float fm, float dt, float amp) {
  for (int it = 0; it < nt; it++) {
    float tmp = M_PI * fm * (it * dt - 1.0 / fm);
    tmp *= tmp;
    wlt[it] = amp * (1.0 - 2.0 * tmp) * std::exp(-tmp);
  }
}

std::vector<float> rickerWavelet(int nt, float fm, float dt, float amp) {
  std::vector<float> wlt(nt);
  for (int it = 0; it < nt; it++) {
    float tmp = M_PI * fm * (it * dt - 1.0 / fm);
    tmp *= tmp;
    wlt[it] = amp * (1.0 - 2.0 * tmp) * std::exp(-tmp);
  }
  return wlt;
}

void rickerWaveletPFwi(float *wlt, int nt, float fm, float dt, float amp) {
  for (int it = 0; it < nt; it++) {
    float tmp = M_PI * fm * ((it - nt / 2) * dt);
    tmp *= tmp;
    wlt[it] = amp * (1.0 - 2.0 * tmp) * std::exp(-tmp);
  }
}

