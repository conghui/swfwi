/*
 * fd4t10s-zjh.cpp
 *
 *  Created on: Mar 5, 2016
 *      Author: rice
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>
#include <cmath>
#include <cstring>
#include <unistd.h>

#define NUM_THREADS 3
#define D 6
const int nt = 2;
const int nx = 20;
const int nz = 20;
const int max_stencil_len = 2 * D + 1;
const int max_strip_4_ldm = 15; /// must be >= max_stencil_len
const int update_strip_len = max_strip_4_ldm - max_stencil_len + 1;

pthread_barrier_t barr;

struct args_t {
  float *p0;
  float *p1;
  float *v;
  int nx;
  int nz;
  int nthreads;
  int tid;
};

float myrand() {
  return static_cast<float>(rand() % 10);
}

void print(const float *A, int nx, int nz) {
  for (int iz = 0; iz < nz ; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      std::printf("%6.1f", A[ix * nz + iz]);
    }
    std::printf("\n");
  }
    std::printf("\n");
}

void print(float **A, int nx, int nz) {
  for (int iz = 0; iz < nz ; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      std::printf("%6.1f", A[ix][iz]);
    }
    std::printf("\n");
  }
    std::printf("\n");
}


void fd4t10s_zjh_2d_vtrans(float *prev_wave, const float *curr_wave, float *next_wave, const float *vel, float *u2, int nx, int nz) {
  float a[6];

  const int d = 6;
  int ix, iz;

  /// Zhang, Jinhai's method
  a[0] = +1.53400796;
  a[1] = +1.78858721;
  a[2] = -0.31660756;
  a[3] = +0.07612173;
  a[4] = -0.01626042;
  a[5] = +0.00216736;

  for (ix = d - 1; ix < nx - (d - 1); ix ++) {
    for (iz = d - 1; iz < nz - (d - 1); iz ++) {
      int curPos = ix * nz + iz;
      u2[curPos] = -4.0 * a[0] * curr_wave[curPos] +
                   a[1] * (curr_wave[curPos - 1]  +  curr_wave[curPos + 1]  +
                           curr_wave[curPos - nz]  +  curr_wave[curPos + nz])  +
                   a[2] * (curr_wave[curPos - 2]  +  curr_wave[curPos + 2]  +
                           curr_wave[curPos - 2 * nz]  +  curr_wave[curPos + 2 * nz])  +
                   a[3] * (curr_wave[curPos - 3]  +  curr_wave[curPos + 3]  +
                           curr_wave[curPos - 3 * nz]  +  curr_wave[curPos + 3 * nz])  +
                   a[4] * (curr_wave[curPos - 4]  +  curr_wave[curPos + 4]  +
                           curr_wave[curPos - 4 * nz]  +  curr_wave[curPos + 4 * nz])  +
                   a[5] * (curr_wave[curPos - 5]  +  curr_wave[curPos + 5]  +
                           curr_wave[curPos - 5 * nz]  +  curr_wave[curPos + 5 * nz]);

    }
  }

  for (ix = d; ix < nx - d; ix++) { /// the range of ix is different from that in previous for loop
    for (iz = d; iz < nz - d; iz++) { /// be careful of the range of iz
      int curPos = ix * nz + iz;
      float curvel = vel[curPos];

      next_wave[curPos] = 2. * curr_wave[curPos] - prev_wave[curPos]  +
                          (1.0f / curvel) * u2[curPos] + /// 2nd order
                          1.0f / 12 * (1.0f / curvel) * (1.0f / curvel) *
                          (u2[curPos - 1] + u2[curPos + 1] + u2[curPos - nz] + u2[curPos + nz] - 4 * u2[curPos]); /// 4th order
    }
  }

}

void fd4t10s_zjh_2d_vtrans_2(float **prev_wave, float **curr_wave, float **next_wave, float **vel, float **u2, int nx, int nz) {
  float a[6];

  const int d = 6;
  int ix, iz;

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
  float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
  std::printf("it: %d, sum of p0: %.2f\n", -1, s);

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
    fd4t10s_zjh_2d_vtrans_2(p0_ptr, p1_ptr, p0_ptr, v_ptr, u2_ptr, nx, nz);
    float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
    sum.push_back(s);
    std::printf("it: %d, sum of p0: %.2f\n", it, s);

    std::swap(p0, p1);
  }

  return sum;
}

void *stencil_wrapper2(void *_arg) {
  struct args_t arg = *(struct args_t *)_arg;
  int nstrip_per_thread = std::ceil(1.0 * arg.nx / arg.nthreads); // # of rows for each threads
  int strip_begin = std::max(nstrip_per_thread * arg.tid - D, 0); // x starts for each threads
  int strip_end =  std::min(nstrip_per_thread * (arg.tid + 1) + D, arg.nx); /// x ends for each threads
  int strip_len = std::max(strip_end - strip_begin, 0);

  int host_istrip;
  if (arg.tid == 0) {
    //printf("original:\n");

    //for (int ix = 0; ix < max_strip_4_ldm; ix++) {
      //for (int iz = 0; iz < nz; iz++) {
        //printf("%4d", arg.p0[ix * nz + iz]);
      //}
      //printf("\n");
    //}
    //printf("\n");
    //printf("copied:\n");
  }

  float ldm_p0[max_strip_4_ldm][nz] = {0};
  float ldm_p1[max_strip_4_ldm][nz] = {0};
  float ldm_p2[max_strip_4_ldm][nz] = {0};
  float ldm_v[max_strip_4_ldm][nz]  = {0};
  float ldm_u2[max_strip_4_ldm][nz] = {0};

  float *ldm_p0_ptr[max_strip_4_ldm];
  float *ldm_p1_ptr[max_strip_4_ldm];
  float *ldm_p2_ptr[max_strip_4_ldm];
  float *ldm_v_ptr[max_strip_4_ldm];
  float *ldm_u2_ptr[max_strip_4_ldm];

  for (host_istrip = 0; host_istrip < max_strip_4_ldm; host_istrip++) {
    ldm_p0_ptr[host_istrip] = ldm_p0[host_istrip];
    ldm_p1_ptr[host_istrip] = ldm_p1[host_istrip];
    ldm_p2_ptr[host_istrip] = ldm_p2[host_istrip];
    ldm_v_ptr[host_istrip] = ldm_v[host_istrip];
    ldm_u2_ptr[host_istrip] = ldm_u2[host_istrip];
  }

  for (host_istrip = 0; host_istrip < std::min(max_strip_4_ldm - update_strip_len, strip_len); host_istrip++) {
    std::memcpy(ldm_p0[host_istrip], arg.p0 + (host_istrip + strip_begin) * nz, nz * sizeof *arg.p0);
    std::memcpy(ldm_p2[host_istrip], arg.p0 + (host_istrip + strip_begin) * nz, nz * sizeof *arg.p0);
    std::memcpy(ldm_p1[host_istrip], arg.p1 + (host_istrip + strip_begin) * nz, nz * sizeof *arg.p1);
    std::memcpy(ldm_v[host_istrip], arg.v + (host_istrip + strip_begin ) * nz, nz * sizeof *arg.v);

    //if (arg.tid == 0 ) {
      //for (int i = 0; i < nz; i++) {
        //printf("%4d", ldm_p0[host_istrip][i]);
      //}
      //printf("\n");
    //}
  }
  pthread_barrier_wait(&barr);

  //if (arg.tid == 0) {
    //printf("nstrip_per_thread: %d\n", nstrip_per_thread);
    //printf("strip_len: %d\n", strip_len);

    ////for (int iz = 0; iz < nz; iz++) {
      ////for (int ix = 0; ix < max_strip_4_ldm; ix++) {
        ////printf("%4d", ldm_p0[ix][iz]);
      ////}
      ////printf("\n");
    ////}
    //print(ldm_p0[0], max_strip_4_ldm, nz);
  //}
  //pthread_barrier_wait(&barr);

  int dstPos = strip_begin + D;

  int watched_thread = 100;
  while (host_istrip < strip_len) {
    int num_strip_to_copy = std::min(update_strip_len, strip_len - host_istrip);
    if (arg.tid == watched_thread) {
      printf("before copy\n");
      print(ldm_p0_ptr, max_strip_4_ldm, nz);
    }

    /// copy the strip
    for (int i = 0; i < num_strip_to_copy; i++) {
      std::memcpy(ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i], arg.p0 + (i + host_istrip + strip_begin) * nz, nz * sizeof *arg.p0);
      std::memcpy(ldm_p2_ptr[max_strip_4_ldm - update_strip_len + i], arg.p0 + (i + host_istrip + strip_begin) * nz, nz * sizeof *arg.p0);
      std::memcpy(ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i], arg.p1 + (i + host_istrip + strip_begin) * nz, nz * sizeof *arg.p1);
      std::memcpy(ldm_v_ptr[max_strip_4_ldm - update_strip_len + i], arg.v +   (i + host_istrip + strip_begin ) * nz, nz * sizeof *arg.v);
    }

    if (arg.tid == watched_thread) {
      printf("before stencil\n");
      print(ldm_p0_ptr, max_strip_4_ldm, nz);
    }

//    usleep(0.1 * 1e6);
//  pthread_barrier_wait(&barr);
//    stencil_kernel_3(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);
    fd4t10s_zjh_2d_vtrans_2(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, ldm_u2_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);

    if (arg.tid == watched_thread) {
      printf("after stencil: %d\n", host_istrip);
      print(ldm_p2_ptr, max_strip_4_ldm, nz);
    }

    //if (host_istrip == 5) {
      //pthread_barrier_wait(&barr);
    //}

    ///copy the data out
    for (int i = 0; i < num_strip_to_copy; i++) {
//      int dst_istrp = host_istrip - D + i;
//      int dst_istrp = host_istrip - update_strip_len + i;
//      int dst_istrp = host_istrip - num_strip_to_copy + i;
      int dst_istrp = dstPos;
      printf("dstpos: %d\n", dstPos);
//      int src_istrp = max_strip_4_ldm - (D + 1) + i;
      int src_istrp = D + i;
      std::memcpy(arg.p0 + (dst_istrp) * nz, ldm_p2_ptr[src_istrp], nz * sizeof *arg.p0);

      dstPos++;
    }

//    usleep(0.1 * 1e6);
    if (arg.tid == watched_thread) {
      printf("global p:\n");
      print(arg.p0, nx, nz);
    }

    /// copy next strip in
    float *p0tmp[update_strip_len];
    float *p1tmp[update_strip_len];
    float *p2tmp[update_strip_len];
    float *vtmp[update_strip_len];
    for (int i = 0; i < update_strip_len; i++) {
      p0tmp[i] = ldm_p0_ptr[i];
      p1tmp[i] = ldm_p1_ptr[i];
      p2tmp[i] = ldm_p2_ptr[i];
      vtmp[i] = ldm_v_ptr[i];
    }
    for (int i = 0; i < max_strip_4_ldm - update_strip_len; i++) {
      ldm_p0_ptr[i] = ldm_p0_ptr[i + update_strip_len];
      ldm_p1_ptr[i] = ldm_p1_ptr[i + update_strip_len];
      ldm_p2_ptr[i] = ldm_p2_ptr[i + update_strip_len];
      ldm_v_ptr[i]  = ldm_v_ptr[i+ update_strip_len];
    }

//    for (int i = max_strip_4_ldm - update_strip_len; i < max_strip_4_ldm; i++) {
    for (int i = 0; i < update_strip_len; i++) {
      ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i] = p0tmp[i];
      ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i] = p1tmp[i];
      ldm_p2_ptr[max_strip_4_ldm - update_strip_len + i] = p2tmp[i];
      ldm_v_ptr[max_strip_4_ldm - update_strip_len + i]  = vtmp[i];
    }


//    int *p0tmp = ldm_p0_ptr[0];
//    int *p1tmp = ldm_p1_ptr[0];
//    int *p2tmp = ldm_p2_ptr[0];
//    int *vtmp = ldm_v_ptr[0];
//    for (int i = 0; i < max_strip_4_ldm - 1; i++) {
//      ldm_p0_ptr[i] = ldm_p0_ptr[i + 1];
//      ldm_p1_ptr[i] = ldm_p1_ptr[i + 1];
//      ldm_p2_ptr[i] = ldm_p2_ptr[i + 1];
//      ldm_v_ptr[i]  = ldm_v_ptr[i+ 1];
//    }
//    ldm_p0_ptr[max_strip_4_ldm - 1] = p0tmp;
//    ldm_p1_ptr[max_strip_4_ldm - 1] = p1tmp;
//    ldm_p2_ptr[max_strip_4_ldm - 1] = p2tmp;
//    ldm_v_ptr[max_strip_4_ldm - 1]  = vtmp;

    host_istrip += num_strip_to_copy;

  }

  return NULL;
}

std::vector<float> run_parallel(float *p0, float *p1, float *v, int nx, int nz, int nt) {
  pthread_t threads[NUM_THREADS];
  struct args_t arg[NUM_THREADS];

  std::vector<float> sum;
  float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
  std::printf("it: %d, sum of p0: %.2f\n", -1, s);

  for (int it = 0; it < nt; it++) {

    for (int tid = 0; tid < NUM_THREADS; tid++) {
      arg[tid].p0 = p0;
      arg[tid].p1 = p1;
      arg[tid].v = v;
      arg[tid].nx = nx;
      arg[tid].nz = nz;
      arg[tid].nthreads = NUM_THREADS;
      arg[tid].tid = tid;
     int rc = pthread_create(&threads[tid], NULL, stencil_wrapper2, (void *)&arg[tid]);

     if (rc){
       printf("ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
     }
    }

    for (int tid = 0; tid < NUM_THREADS; tid++) {
      pthread_join(threads[tid], NULL);
    }
    float s = std::accumulate(p0, p0 + nx * nz, 0.0f);
    sum.push_back(s);
    std::printf("it: %d, sum of p0: %.2f\n", it, s);
    std::swap(p0, p1);
  }

  return sum;
}

int main(int argc, char *argv[])
{
  std::vector<float> p0(nx * nz);
  std::vector<float> p1(nx * nz);
  std::vector<float> v(nx * nz);

  if(pthread_barrier_init(&barr, NULL, NUM_THREADS))
  {
    printf("Could not create a barrier\n");
    return -1;
  }

  std::generate(p0.begin(), p0.end(), myrand);
  std::generate(p1.begin(), p1.end(), myrand);
  std::fill(v.begin(), v.end(), 1500);

  std::printf("p0:\n");
  print(&p0[0], nx, nz);

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
    print(&tmp_p0[0], nx, nz);
  }

  std::vector<float> sum_parallel;
  {
    std::vector<float> tmp_p0 = p0;
    std::vector<float> tmp_p1 = p1;
    std::vector<float> tmp_v = v;
    sum_parallel = run_parallel(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
    print(&tmp_p0[0], nx, nz);
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

