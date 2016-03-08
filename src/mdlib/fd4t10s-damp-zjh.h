/*
 * fd4t10s-damp-zjh.h
 *
 *  Created on: Mar 4, 2016
 *      Author: rice
 */

#ifndef SRC_MDLIB_FD4T10S_DAMP_ZJH_H_
#define SRC_MDLIB_FD4T10S_DAMP_ZJH_H_

/**
 * please note that the velocity is not transformed
 */
void fd4t10s_damp_zjh_2d(float *prev_wave, const float *curr_wave, const float *vel, int nx, int nz, int nb, float dx, float dt);

#endif /* SRC_MDLIB_FD4T10S_DAMP_ZJH_H_ */
