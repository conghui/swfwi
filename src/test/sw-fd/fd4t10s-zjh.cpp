/*
 * fd4t10s-zjh.cpp
 *
 *  Created on: Mar 5, 2016
 *      Author: rice
 */

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "fd-arg.h"

extern "C" {
#include "athread.h"
void slave_func(void *);
}

static inline unsigned long rpcc()
{
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}

float myrand() {
  return static_cast<float>(rand() % 10);
}

void print(const float *A, int nx, int nz) {
  for (int iz = 0; iz < nz ; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      std::printf("%5.0f", A[ix * nz + iz]);
    }
    std::printf("\n");
  }
    std::printf("\n");
}

void print2(float **A, int nx, int nz) {
  for (int iz = 0; iz < nz ; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      std::printf("%6.1f", A[ix][iz]);
    }
    std::printf("\n");
  }
    std::printf("\n");
}

void fd4t10s_zjh_2d_vtrans_2(float **prev_wave, float **curr_wave, float **next_wave, float **vel, float **u2, int nx, int nz) {
  float a[6];
  int ix, iz;
  const int d = 6;

  /// Zhang, Jinhai's method
  a[0] = +1.53400796;
  a[1] = +1.78858721;
  a[2] = -0.31660756;
  a[3] = +0.07612173;
  a[4] = -0.01626042;
  a[5] = +0.00216736;

  for (ix = d - 1; ix < nx - (d - 1); ix ++) {
    for (iz = d - 1; iz < nz - (d - 1); iz ++) {
      u2[ix][iz] = -4.0 * a[0] * curr_wave[ix][iz] +
                   a[1] * (curr_wave[ix][iz-1] + curr_wave[ix][iz+1] + curr_wave[ix-1][iz] + curr_wave[ix+1][iz]) +
                   a[2] * (curr_wave[ix][iz-2] + curr_wave[ix][iz+2] + curr_wave[ix-2][iz] + curr_wave[ix+2][iz]) +
                   a[3] * (curr_wave[ix][iz-3] + curr_wave[ix][iz+3] + curr_wave[ix-3][iz] + curr_wave[ix+3][iz]) +
                   a[4] * (curr_wave[ix][iz-4] + curr_wave[ix][iz+4] + curr_wave[ix-4][iz] + curr_wave[ix+4][iz]) +
                   a[5] * (curr_wave[ix][iz-5] + curr_wave[ix][iz+5] + curr_wave[ix-5][iz] + curr_wave[ix+5][iz]);
    }
  }

  for (ix = d; ix < nx - d; ix++) { /// the range of ix is different from that in previous for loop
    for (iz = d; iz < nz - d; iz++) { /// be careful of the range of iz
      float curvel = vel[ix][iz];
      next_wave[ix][iz] = 2. * curr_wave[ix][iz] - prev_wave[ix][iz]  +
                          (1.0f / curvel) * u2[ix][iz] + /// 2nd order
                          1.0f / 12 * (1.0f / curvel) * (1.0f / curvel) *
                          (u2[ix][iz-1] + u2[ix][iz+1] + u2[ix-1][iz] + u2[ix+1][iz] - 4 * u2[ix][iz]); /// 4th order
    }
  }
}

std::vector<float> run_serial(float *p0, float *p1, float *v, int nx, int nz, int nt) {
  float *p0_ptr[nx];
  float *p1_ptr[nx];
  float *v_ptr[nx];
  float *u2_ptr[nx];

  std::vector<float> sum;
  std::vector<float> u2(nx*nz, 0);
  //float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
  //std::printf("it: %d, sum of p0: %.2f\n", -1, s);

  for (int it = 0; it < nt; it++) {
    for (int ix = 0; ix < nx; ix++) {
      p0_ptr[ix] = &p0[ix * nz];
      p1_ptr[ix] = &p1[ix * nz];
      v_ptr[ix] = &v[ix * nz];
      u2_ptr[ix] = &u2[ix*nz];
    }
    //stencil_kernel(p0, p1, v, nx, nz);
//    stencil_kernel_2(p0_ptr, p1_ptr, v_ptr, nx, nz);
//    print(p0_ptr[0], nx, nz);
//    fd4t10s_zjh_2d_vtrans(p0, p1, p0, v, &u2[0], nx, nz);
    int st = rpcc();
    fd4t10s_zjh_2d_vtrans_2(p0_ptr, p1_ptr, p0_ptr, v_ptr, u2_ptr, nx, nz);
    int ed = rpcc();

    float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
    sum.push_back(s);
    std::printf("it: %2d, sum of p0: %.2f, count: %d\n", it, s, ed - st);

    std::swap(p0, p1);
  }

  return sum;
}

std::vector<float> run_parallel(float *p0, float *p1, float *v, int nx, int nz, int nt) {
  struct args_t arg[NUM_THREADS];

  std::vector<float> sum;
  //float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
  //std::printf("it: %d, sum of p0: %.2f\n", -1, s);

  athread_init();

  //struct args_t arg;
  //arg.p0 = p0;

  for (int it = 0; it < nt; it++) {

    int st = rpcc();
    for (int tid = 0; tid < NUM_THREADS; tid++) {
      arg[tid].p0 = p0;
      arg[tid].p1 = p1;
      arg[tid].v = v;
      arg[tid].nx = nx;
      arg[tid].nz = nz;
      arg[tid].nthreads = NUM_THREADS;
     int rc = __real_athread_create(tid, (void *)slave_func, &arg[tid]);

     //if (rc != 0){
       //printf("ERROR; return code from pthread_create() is %d\n", rc);
       //exit(-1);
     //}
    }

    for (int tid = 0; tid < NUM_THREADS; tid++) {
      //pthread_join(threads[tid], NULL);
      if (athread_wait(tid) != 0) {
        printf("error\n");
        exit(0);
      }
    }
    for (int tid = 0; tid < NUM_THREADS; tid++) {
      //pthread_join(threads[tid], NULL);
      if (athread_end(tid) != 0) {
        printf("error\n");
        exit(0);
      }
    }
    int ed = rpcc();
    //__real_athread_spawn((void *)slave_func, &arg);

    float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
    sum.push_back(s);
    std::printf("it: %2d, sum of p0: %.2f, count: %d\n", it, s, ed - st);
    std::swap(p0, p1);
  }

  return sum;
}

int main(int argc, char *argv[])
{
  int nx = NX;
  int nz = NZ;
  int nt = NT;
  std::vector<float> p0(nx * nz);
  std::vector<float> p1(nx * nz);
  std::vector<float> v(nx * nz);

  std::generate(p0.begin(), p0.end(), myrand);
  std::generate(p1.begin(), p1.end(), myrand);
  std::fill(v.begin(), v.end(), 1500);

  std::printf("p0:\n");
  //print(&p0[0], nx, nz);

  //std::printf("\np1:\n");
  //print(&p1[0], nx, nz);

  //std::printf("\nv:\n");
  //print(&v[0], nx, nz);

  std::vector<float> sum_serial;
  {
    std::vector<float> tmp_p0 = p0;
    std::vector<float> tmp_p1 = p1;
    std::vector<float> tmp_v = v;
    sum_serial = run_serial(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
    //print(&tmp_p0[0], nx, nz);
  }

  std::vector<float> sum_parallel;
  {
    std::vector<float> tmp_p0 = p0;
    std::vector<float> tmp_p1 = p1;
    std::vector<float> tmp_v = v;
    sum_parallel = run_parallel(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
    //print(&tmp_p0[0], nx, nz);
  }

  for (int i = 0; i < nt; i++) {
    if (sum_serial[i] != sum_parallel[i]) {
      printf("error\n");
      exit(0);
    }
  }

  printf("test pass!\n");

  return 0;
}

