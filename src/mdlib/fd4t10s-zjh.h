/*
 * fd4t10s-zjh.h
 *
 *  Created on: Mar 5, 2016
 *      Author: rice
 */

#ifndef SRC_MDLIB_FD4T10S_ZJH_H_
#define SRC_MDLIB_FD4T10S_ZJH_H_

void fd4t10s_zjh_2d_vtrans(float *prev_wave, const float *curr_wave, const float *vel, int nx, int nz);
void fd4t10s_zjh_2d_vtrans_lap_illum(float *lap, float *illum, float *prev_wave, const float *curr_wave, const float *vel, int nx, int nz);

void fd4t10s_zjh_2d_vnotrans(float *prev_wave, const float *curr_wave, const float *vel, int nx, int nz, float dx, float dt);
void fd4t10s_zjh_2d_vnotrans_lap_illum(float *lap, float *illum, float *prev_wave, const float *curr_wave, const float *vel, int nx, int nz, float dx, float dt);

#endif /* SRC_MDLIB_FD4T10S_ZJH_H_ */
