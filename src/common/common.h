/*
 * common.h
 *
 *  Created on: Feb 24, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_COMMON_H_
#define SRC_COMMON_COMMON_H_

#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

template <typename T>
void vectorMinus(const std::vector<T> &dobs, const std::vector<T> &dcal, std::vector<T> &vsrc) {
  std::transform(dobs.begin(), dobs.end(), dcal.begin(), vsrc.begin(), std::minus<T>());
}

template <typename E>
float velRecover(float vel, float dx, float dt) {
  float t = dx * dx / (dt * dt * vel);
  return std::sqrt(t);
}

template <typename E>
float velTrans(float vel, float dx, float dt) {
  return dx * dx / (dt * dt * vel * vel);
}

template <typename T>
bool abs_less(T a, T b) {
  return std::abs(a) < std::abs(b);
}

template <typename T>
float variance(const T *A_Begin, const T *A_End, const T *B_Begin) {
  float v = 0;

  for (; A_Begin != A_End; ++A_Begin, ++B_Begin) {
    v += (*A_Begin - *B_Begin) * (*A_Begin - *B_Begin);
  }

  return v;
}


void matrix_transpose(float *matrix, float *trans, int n1, int n2);
void step_forward(const float *p0, const float *p1, float *p2, const float *vv, float dtz, float dtx, int nz, int nx);
void step_backward(float *illum, float *lap, const float *p0, const float *p1, float *p2, const float *vv, float dtz, float dtx, int nz, int nx);

void add_source(float *p, const float *source, const int *sxz, int ns, int nz, int nb, bool add);

void record_seis(float *seis_it, const int *gxz, const float *p, int ng, int nz, int nb);
void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz);

void rw_bndr(float *bndr, float *p, int nz, int nx, bool write);

void cal_residuals(const float *dcal, const float *dobs, float *dres, int ng);


void cal_gradient(float *grad, const float *lap, const float *gp, int nz, int nx);


void scale_gradient(float *grad, const float *vv, const float *illum, int nz, int nx, bool precon);

float cal_objective(float *dres, int ng);

float cal_beta(const float *g0, const float *g1, const float *cg, int nz, int nx);

void cal_conjgrad(const float *g1, float *cg, float beta, int nz, int nx);


float cal_epsilon(const float *vv, const float *cg, int nz, int nx);

void cal_vtmp(float *vtmp, const float *vv, const float *cg, float epsil, int nz, int nx);


void sum_alpha12(float *alpha1, float *alpha2, const float *dcaltmp, const float *dobs, const float *derr, int ng);

void bell_smoothz(const float *g, float *smg, int rbell, int nz, int nx);
void bell_smoothx(const float *g, float *smg, int rbell, int nz, int nx);
float cal_alpha(const float *alpha1, const float *alpha2, float epsil, int ng);
void update_vel(float *vv, const float *cg, float alpha, int nz, int nx);

#endif /* SRC_COMMON_COMMON_H_ */
