/*
 * ricker-wavelet.h
 *
 *  Created on: Feb 24, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_RICKER_WAVELET_H_
#define SRC_COMMON_RICKER_WAVELET_H_

#include <vector>

void rickerWavelet(float *wlt, int nt, float fm, float dt, float amp);
void rickerWaveletPFwi(float *wlt, int nt, float fm, float dt, float amp);
std::vector<float> rickerWavelet(int nt, float fm, float dt, float amp);

#endif /* SRC_COMMON_RICKER_WAVELET_H_ */
