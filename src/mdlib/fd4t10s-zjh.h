/*
 * fd4t10s-zjh.h
 *
 *  Created on: Mar 5, 2016
 *      Author: rice
 */

#ifndef SRC_MDLIB_FD4T10S_ZJH_H_
#define SRC_MDLIB_FD4T10S_ZJH_H_

void fd4t10s_zjh_2d(float *prev_wave, const float *curr_wave, const float *vel, int nx, int nz, int nb, float dx, float dt);

#endif /* SRC_MDLIB_FD4T10S_ZJH_H_ */
