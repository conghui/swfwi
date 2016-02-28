/*
 * common.cpp
 *
 *  Created on: Feb 24, 2016
 *      Author: rice
 */

extern "C" {
#include <rsf.h>
}
#include "common.h"

void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{

#pragma omp parallel for
  for (int i2 = 0; i2 < n2; i2++) {
    for (int i1 = 0; i1 < n1; i1++) {
      trans[i2 + n2 * i1] = matrix[i1 + n1 * i2];
    }
  }
}


void step_forward(const float *p0, const float *p1, float *p2, const float *vv, float dtz, float dtx, int nz, int nx)
/*< forward modeling step, Clayton-Enquist ABC incorporated >*/
{

#pragma omp parallel for
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      float v1 = vv[ix * nz + iz] * dtz;
      v1 = v1 * v1;
      float v2 = vv[ix * nz + iz] * dtx;
      v2 = v2 * v2;
      float diff1 = -2.0 * p1[ix * nz + iz];
      float diff2 = diff1;

      diff1 += (iz - 1 >= 0) ? p1[ix * nz + (iz - 1)] : 0.0;
      diff1 += (iz + 1 < nz) ? p1[ix * nz + (iz + 1)] : 0.0;
      diff2 += (ix - 1 >= 0) ? p1[(ix - 1) * nz + iz] : 0.0;
      diff2 += (ix + 1 < nx) ? p1[(ix + 1) * nz + iz] : 0.0;
      diff1 *= v1;
      diff2 *= v2;
      p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
    }
  }

#pragma omp parallel for
  for (int ix = 1; ix < nx - 1; ix++) {
    /* top boundary */
    /*
    iz=0;
    diff1= (p1[ix][iz+1]-p1[ix][iz])-
    (p0[ix][iz+1]-p0[ix][iz]);
    diff2= c21*(p1[ix-1][iz]+p1[ix+1][iz]) +
    c22*(p1[ix-2][iz]+p1[ix+2][iz]) +
    c20*p1[ix][iz];
    diff1*=sqrtf(vv[ix][iz])/dz;
    diff2*=vv[ix][iz]/(2.0*dx*dx);
    p2[ix][iz]=2*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
     */
    /* bottom boundary */
    int iz = nz - 1;
    float v1 = vv[ix * nz + iz] * dtz;
    float v2 = vv[ix * nz + iz] * dtx;
    float diff1 = -(p1[ix * nz + iz] - p1[ix * nz + (iz - 1)]) + (p0[ix * nz + iz] - p0[ix * nz + iz - 1]);
    float diff2 = p1[(ix - 1) * nz + iz] - 2.0 * p1[ix * nz + iz] + p1[(ix + 1) * nz + iz];
    diff1 *= v1;
    diff2 *= 0.5 * v2 * v2;
    p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
  }

#pragma omp parallel for
  for (int iz = 1; iz < nz - 1; iz++) {
    /* left boundary */
    int ix = 0;
    float v1 = vv[ix * nz + iz] * dtz;
    float v2 = vv[ix * nz + iz] * dtx;
    float diff1 = p1[ix * nz + (iz - 1)] - 2.0 * p1[ix * nz + iz] + p1[ix * nz + (iz + 1)];
    float diff2 = (p1[(ix + 1) * nz + iz] - p1[ix * nz + iz]) - (p0[(ix + 1) * nz + iz] - p0[ix * nz + iz]);
    diff1 *= 0.5 * v1 * v1;
    diff2 *= v2;
    p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;

    /* right boundary */
    ix = nx - 1;
    v1 = vv[ix * nz + iz] * dtz;
    v2 = vv[ix * nz + iz] * dtx;
    diff1 = p1[ix * nz + (iz - 1)] - 2.0 * p1[ix * nz + iz] + p1[ix * nz + (iz + 1)];
    diff2 = -(p1[ix * nz + iz] - p1[(ix - 1) * nz + iz]) + (p0[ix * nz + iz] - p0[(ix - 1) * nz + iz]);
    diff1 *= 0.5 * v1 * v1;
    diff2 *= v2;
    p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
  }
}

void step_backward(float *illum, float *lap, const float *p0, const float *p1, float *p2, const float *vv, float dtz, float dtx, int nz, int nx)
/*< step backward >*/
{
#pragma omp parallel for
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      float v1 = vv[ix * nz + iz] * dtz;
      v1 = v1 * v1;
      float v2 = vv[ix * nz + iz] * dtx;
      v2 = v2 * v2;
      float diff1 = -2.0 * p1[ix * nz + iz];
      float diff2 = diff1;
      diff1 += (iz - 1 >= 0) ? p1[ix * nz + (iz - 1)] : 0.0;
      diff1 += (iz + 1 < nz) ? p1[ix * nz + (iz + 1)] : 0.0;
      diff2 += (ix - 1 >= 0) ? p1[(ix - 1) * nz + iz] : 0.0;
      diff2 += (ix + 1 < nx) ? p1[(ix + 1) * nz + iz] : 0.0;
      lap[ix * nz + iz] = diff1 + diff2;
      diff1 *= v1;
      diff2 *= v2;
      p2[ix * nz + iz] = 2.0 * p1[ix * nz + iz] - p0[ix * nz + iz] + diff1 + diff2;
      illum[ix * nz + iz] += p1[ix * nz + iz] * p1[ix * nz + iz];
    }
  }
}

void add_source(float *p, const float *source, const int *sxz, int ns, int nz, int nb, bool add)
/*< add/subtract seismic sources >*/
{
  if (add) {
    for (int is = 0; is < ns; is++) {
      int sx = sxz[is] / nz + nb;
      int sz = sxz[is] % nz;
      p[sx * (nz + nb) + sz] += source[is];
    }
  } else {
    for (int is = 0; is < ns; is++) {
      int sx = sxz[is] / nz + nb;
      int sz = sxz[is] % nz;
      p[sx * (nz + nb) + sz] -= source[is];
    }
  }
}

void record_seis(float *seis_it, const int *gxz, const float *p, int ng, int nz, int nb)
/*< record seismogram at time it into a vector length of ng >*/
{
#pragma omp parallel for
  for (int ig = 0; ig < ng; ig++) {
    int gx = gxz[ig] / nz + nb;
    int gz = gxz[ig] % nz;
    seis_it[ig] = p[gx * (nz + nb) + gz];
  }
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz)
/*< shot/geophone position initialize >*/
{
  for (int is = 0; is < ns; is++) {
    int sz = szbeg + is * jsz;
    int sx = sxbeg + is * jsx;
    sxz[is] = sz + nz * sx;
  }
}

void rw_bndr(float *bndr, float *p, int nz, int nx, bool write)
/*< if write==true, write/save boundaries out of variables;
 else  read boundaries into variables (for 2nd order FD) >*/
{
  if (write) {
    for (int i = 0; i < nz; i++) {
      bndr[i] = p[i];
      bndr[i + nz] = p[(nx - 1) * nz + i];
    }
    for (int i = 0; i < nx; i++) {
      bndr[i + 2 * nz] = p[i * nz + (nz - 1)];
    }
  } else {
    for (int i = 0; i < nz; i++) {
      p[i] = bndr[i];
      p[(nx - 1) * nz + i] = bndr[i + nz];
    }
    for (int i = 0; i < nx; i++) {
      p[i * nz + (nz - 1)] = bndr[i + 2 * nz];
    }
  }
}

void cal_residuals(const float *dcal, const float *dobs, float *dres, int ng)
/*< calculate residual >*/
{
  for (int ig = 0; ig < ng; ig++) {
    dres[ig] = dcal[ig] - dobs[ig];
  }
}

void cal_gradient(float *grad, const float *lap, const float *gp, int nz, int nx)
/*< calculate gradient >*/
{
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      int idx = ix * nz + iz;
      grad[idx] += lap[idx] * gp[idx];
    }
  }
}

void scale_gradient(float *grad, const float *vv, const float *illum, int nz, int nx, bool precon)
/*< scale gradient >*/
{
  for (int ix = 1; ix < nx - 1; ix++) {
    for (int iz = 1; iz < nz - 1; iz++) {
      float a = vv[ix * nz + iz];
      if (precon) {
        a *= sqrtf(illum[ix * nz + iz] + SF_EPS);  /*precondition with residual wavefield illumination*/
      }
      grad[ix * nz + iz] *= 2.0 / a;
    }
  }

  for (int ix = 0; ix < nx; ix++) {
    grad[ix * nz + 0] = grad[ix * nz + 1];
    grad[ix * nz + nz - 1] = grad[ix * nz + nz - 2];
  }

  for (int iz = 0; iz < nz; iz++) {
    grad[iz] = grad[1 * nz + iz];
    grad[(nx - 1) * nz + iz] = grad[(nx - 2) * nz + iz];
  }
}

float cal_objective(float *dres, int ng)
/*< calculate the value of objective function >*/
{
  float obj = 0;

  for (int i = 0; i < ng; i++) {
    float a = dres[i];
    obj += a * a;
  }
  return obj;
}

float cal_beta(const float *g0, const float *g1, const float *cg, int nz, int nx)
/*< calculate beta >*/
{
  float a = 0;
  float b = 0;
  float c = 0;
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      int idx = ix * nz + iz;
      a += g1[idx] * (g1[idx] - g0[idx]); // numerator of HS
      b += cg[idx] * (g1[idx] - g0[idx]); // denominator of HS,DY
      c += g1[idx] * g1[idx]; // numerator of DY
    }
  }

  float beta_HS = (fabsf(b) > 0) ? (a / b) : 0.0;
  float beta_DY = (fabsf(b) > 0) ? (c / b) : 0.0;
  return SF_MAX(0.0, SF_MIN(beta_HS, beta_DY));
}


void cal_conjgrad(const float *g1, float *cg, float beta, int nz, int nx)
/*< calculate conjugate gradient >*/
{

  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      int idx = ix * nz + iz;
      cg[idx] = -g1[idx] + beta * cg[idx];
    }
  }
}

float cal_epsilon(const float *vv, const float *cg, int nz, int nx)
/*< calculate epsilcon >*/
{
  float vvmax = 0.0;
  float cgmax = 0.0;

  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      int idx = ix * nz + iz;
      vvmax = SF_MAX(vvmax, fabsf(vv[idx]));
      cgmax = SF_MAX(cgmax, fabsf(cg[idx]));
    }
  }

  return 0.01 * vvmax / (cgmax + SF_EPS);
}

void cal_vtmp(float *vtmp, const float *vv, const float *cg, float epsil, int nz, int nx)
/*< calculate temporary velcity >*/
{

  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      int idx = ix * nz + iz;
      vtmp[idx] = vv[idx] + epsil * cg[idx];
    }
  }
}

void sum_alpha12(float *alpha1, float *alpha2, const float *dcaltmp, const float *dobs, const float *derr, int ng)
/*< calculate numerator and denominator of alpha >*/
{
  for (int ig = 0; ig < ng; ig++) {
    float c = derr[ig];
    float a = dobs[ig] + c; /* since f(mk)-dobs[id]=derr[id], thus f(mk)=b+c; */
    float b = dcaltmp[ig] - a; /* f(mk+epsil*cg)-f(mk) */
    alpha1[ig] -= b * c;
    alpha2[ig] += b * b;
  }
}


float cal_alpha(const float *alpha1, const float *alpha2, float epsil, int ng)
/*< calculate alpha >*/
{

  float a = 0;
  float b = 0;
  for (int ig = 0; ig < ng; ig++) {
    a += alpha1[ig];
    b += alpha2[ig];
  }

  return (a * epsil / (b + SF_EPS));
}

void update_vel(float *vv, const float *cg, float alpha, int nz, int nx)
/*< update velcity >*/
{

  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      int idx = ix * nz + iz;
      vv[idx] += alpha * cg[idx];
    }
  }
}

void bell_smoothz(const float *g, float *smg, int rbell, int nz, int nx)
/*< gaussian bell smoothing for z-axis >*/
{
  for (int ix = 0; ix < nx; ix++)
    for (int iz = 0; iz < nz; iz++) {
      float s = 0.0;
      for (int i = -rbell; i <= rbell; i++) if (iz + i >= 0 && iz + i < nz) {
        s += expf(-(2.0 * i * i) / rbell) * g[ix * nz + iz + i];
      }
      smg[ix * nz + iz] = s;
    }
}

void bell_smoothx(const float *g, float *smg, int rbell, int nz, int nx)
/*< gaussian bell smoothing for x-axis >*/
{
  for (int ix = 0; ix < nx; ix++)
    for (int iz = 0; iz < nz; iz++) {
      float s = 0.0;
      for (int i = -rbell; i <= rbell; i++) if (ix + i >= 0 && ix + i < nx) {
        s += expf(-(2.0 * i * i) / rbell) * g[(ix + i) * nz + iz];
      }
      smg[ix * nz + iz] = s;
    }
}
