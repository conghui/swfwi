#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "slave.h"
#include "fd-arg.h"

#define max_stencil_len  (2 * (D) + 1)
#define max_strip_4_ldm  (16)
#define update_strip_len  ((max_strip_4_ldm) - (max_stencil_len) + 1)

#if update_strip_len < 1
#error update_strip_len must >= 1
#endif

#define MAX(a, b) (a) > (b) ? (a) : (b)
#define MIN(a, b) (a) < (b) ? (a) : (b)


__thread_local float ldm_p0[max_strip_4_ldm + update_strip_len][NZ];
__thread_local float ldm_p1[max_strip_4_ldm + update_strip_len][NZ];
__thread_local float ldm_v[max_strip_4_ldm  + update_strip_len][NZ];
__thread_local float ldm_u2[max_strip_4_ldm][NZ];

__thread_local float *ldm_p0_ptr[max_strip_4_ldm + update_strip_len];
__thread_local float *ldm_p1_ptr[max_strip_4_ldm + update_strip_len];
__thread_local float *ldm_v_ptr[max_strip_4_ldm  + update_strip_len];
__thread_local float *ldm_u2_ptr[max_strip_4_ldm];
__thread_local float *p0tmp[update_strip_len];
__thread_local float *p1tmp[update_strip_len];
__thread_local float *vtmp[update_strip_len];

__thread_local float ldm_p2[update_strip_len][NZ-2*D];
__thread_local float *ldm_p2_ptr[update_strip_len];

__thread_local volatile int get_reply_p0[update_strip_len];
__thread_local volatile int get_reply_p1[update_strip_len];
__thread_local volatile int get_reply_v[update_strip_len];
__thread_local volatile int put_reply[update_strip_len];

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

void fd4t10s_zjh_2d_vtrans_3(float **prev_wave, float **curr_wave, float **next_wave, float **vel, float **u2, int nx, int nz) {
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
      next_wave[ix-d][iz-d] = 2. * curr_wave[ix][iz] - prev_wave[ix][iz]  +
                          (1.0f / curvel) * u2[ix][iz] + /// 2nd order
                          1.0f / 12 * (1.0f / curvel) * (1.0f / curvel) *
                          (u2[ix][iz-1] + u2[ix][iz+1] + u2[ix-1][iz] + u2[ix+1][iz] - 4 * u2[ix][iz]); /// 4th order
    }
  }
}

void print(float **A, int nx, int nz) {
  for (int iz = 0; iz < nz; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      printf("%5.0f", A[ix][iz]);
    }
    printf("\n");
  }
}

void print1d(float *A, int nx, int nz) {
  for (int iz = 0; iz < nz; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      printf("%5.0f", A[ix * nz + iz]);
    }
    printf("\n");
  }
}

void func3(void *_arg) {
  volatile unsigned int get_reply, put_reply;
  struct args_t arg     = *(struct args_t *)_arg;
  int tid               = athread_get_id(-1);
  int nstrip_per_thread = (arg.nx + arg.nthreads - 1) / arg.nthreads; // # of rows for each threads
  int strip_begin       = MAX(nstrip_per_thread * tid - D, 0); // x starts for each threads
  int strip_end         = MIN(nstrip_per_thread * (tid + 1) + D, arg.nx); /// x ends for each threads
  int strip_len         = MAX(strip_end - strip_begin, 0);
  int host_istrip;
  int nz = NZ;

  int debug = tid == 100;
  if (debug) printf("in function: %s\n", __PRETTY_FUNCTION__);

  for (host_istrip = 0; host_istrip < max_strip_4_ldm; host_istrip++) {
    ldm_p0_ptr[host_istrip] = ldm_p0[host_istrip];
    ldm_p1_ptr[host_istrip] = ldm_p1[host_istrip];
    ldm_v_ptr[host_istrip] = ldm_v[host_istrip];
    ldm_u2_ptr[host_istrip] = ldm_u2[host_istrip];
  }

  for (int i = 0; i < update_strip_len; i++) {
    ldm_p2_ptr[i] = ldm_p2[i];
  }

  int num_strip_to_copy = MIN(max_strip_4_ldm - update_strip_len, strip_len);
  get_reply = 0;
  for (host_istrip = 0; host_istrip < num_strip_to_copy; host_istrip++) {
    athread_get(PE_MODE, arg.p0+(host_istrip+strip_begin)*nz, ldm_p0[host_istrip], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.p1+(host_istrip+strip_begin)*nz, ldm_p1[host_istrip], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.v+(host_istrip+strip_begin)*nz,  ldm_v[host_istrip],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
  }
  while (get_reply != 3 * num_strip_to_copy);

  if (debug) printf("numthread:%d\n", arg.nthreads);

  int dstPos = strip_begin + D;
  while (host_istrip < strip_len) {
    num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);

    get_reply = 0;
    for (int i = 0; i < num_strip_to_copy; i++) {
      athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
      athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
      athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm-update_strip_len+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
    }
    while (get_reply != 3 * num_strip_to_copy);

    if (debug) {
      printf("before stencil\n");
      print(ldm_p0_ptr, max_strip_4_ldm, nz);
    }

    fd4t10s_zjh_2d_vtrans_3(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, ldm_u2_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);

    if (debug) {
      printf("after stencil: %d\n", host_istrip);
      print(ldm_p2_ptr, update_strip_len, nz-2*D);
    }


    ///copy the data out
    put_reply = 0;
    for (int i = 0; i < num_strip_to_copy; i++) {
      int dst_istrp = dstPos;
      int src_istrp = i;
      athread_put(PE_MODE, ldm_p2_ptr[src_istrp], arg.p0+(dst_istrp)*nz + D, (nz-2*D)*sizeof *arg.p0, (void *)&put_reply, 0, 0);
      dstPos++;
    }
    while (put_reply != num_strip_to_copy);

    /// copy next strip in
    for (int i = 0; i < update_strip_len; i++) {
      p0tmp[i] = ldm_p0_ptr[i];
      p1tmp[i] = ldm_p1_ptr[i];
      vtmp[i] = ldm_v_ptr[i];
    }
    for (int i = 0; i < max_strip_4_ldm - update_strip_len; i++) {
      ldm_p0_ptr[i] = ldm_p0_ptr[i + update_strip_len];
      ldm_p1_ptr[i] = ldm_p1_ptr[i + update_strip_len];
      ldm_v_ptr[i]  = ldm_v_ptr[i+ update_strip_len];
    }

    for (int i = 0; i < update_strip_len; i++) {
      ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i] = p0tmp[i];
      ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i] = p1tmp[i];
      ldm_v_ptr[max_strip_4_ldm - update_strip_len + i]  = vtmp[i];
    }

    host_istrip += num_strip_to_copy;

  }

  if (debug) {
    printf("tid: %d, exit %s\n", tid, __PRETTY_FUNCTION__);
  }

}

void func(void *_arg) {
  volatile unsigned int get_reply;
  struct args_t arg     = *(struct args_t *)_arg;
  int tid               = athread_get_id(-1);
  int nstrip_per_thread = (arg.nx + arg.nthreads - 1) / arg.nthreads; // # of rows for each threads
  int strip_begin       = MAX(nstrip_per_thread * tid - D, 0); // x starts for each threads
  int strip_end         = MIN(nstrip_per_thread * (tid + 1) + D, arg.nx); /// x ends for each threads
  int strip_len         = MAX(strip_end - strip_begin, 0);
  int host_istrip;
  int nz = NZ;

  int debug = tid == 1000;
  if (debug) printf("in function: %s\n", __PRETTY_FUNCTION__);

  for (int i = 0; i < max_strip_4_ldm + update_strip_len; i++) {
    ldm_p0_ptr[i] = ldm_p0[i];
    ldm_p1_ptr[i] = ldm_p1[i];
    ldm_v_ptr[i] = ldm_v[i];
  }

  for (int i = 0; i < max_strip_4_ldm; i++) {
    ldm_u2_ptr[i] = ldm_u2[i];
  }

  for (int i = 0; i < update_strip_len; i++) {
    ldm_p2_ptr[i] = ldm_p2[i];
  }

  int num_strip_to_copy = MIN(max_strip_4_ldm - update_strip_len, strip_len);
  get_reply = 0;
  for (host_istrip = 0; host_istrip < num_strip_to_copy; host_istrip++) {
    athread_get(PE_MODE, arg.p0+(host_istrip+strip_begin)*nz, ldm_p0[host_istrip], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.p1+(host_istrip+strip_begin)*nz, ldm_p1[host_istrip], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.v+(host_istrip+strip_begin)*nz,  ldm_v[host_istrip],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
  }
  while (get_reply != 3 * num_strip_to_copy);

  if (debug) printf("numthread:%d\n", arg.nthreads);

  num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);

  get_reply = 0;
  for (int i = 0; i < num_strip_to_copy; i++) {
    athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm-update_strip_len+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
  }
  while (get_reply != 3 * num_strip_to_copy);

  host_istrip += num_strip_to_copy;
  int dstPos = strip_begin + D;
  int cur_nx = host_istrip;

  do {
    int next_num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);

    for (int i= 0; i < update_strip_len; i++) {
      get_reply_p0[i] = 0;
      get_reply_p1[i] = 0;
      get_reply_v[i] = 0;
      put_reply[i] = 0;
    }
    for (int i = 0; i < next_num_strip_to_copy; i++) {
      athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm+i], nz*sizeof *arg.p0, (void *)&get_reply_p0[i], 0, 0, 0);
      athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm+i], nz*sizeof *arg.p1, (void *)&get_reply_p1[i], 0, 0, 0);
      athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm+i],  nz*sizeof *arg.v,  (void *)&get_reply_v[i], 0, 0, 0);
    }

    /*while (get_reply != 3 * next_num_strip_to_copy);*/
    fd4t10s_zjh_2d_vtrans_3(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, ldm_u2_ptr, cur_nx, nz);

    ///copy the data out
    for (int i = 0; i < num_strip_to_copy; i++) {
      int dst_istrp = dstPos;
      int src_istrp = i;
      athread_put(PE_MODE, ldm_p2_ptr[src_istrp], arg.p0+(dst_istrp)*nz + D, (nz-2*D)*sizeof *arg.p0, (void *)&put_reply[i], 0, 0);
      dstPos++;
    }

    for (int i = 0; i < update_strip_len; i++) {
      p0tmp[i] = ldm_p0_ptr[i];
      p1tmp[i] = ldm_p1_ptr[i];
      vtmp[i] = ldm_v_ptr[i];
    }
    for (int i = 0; i < max_strip_4_ldm; i++) {
      ldm_p0_ptr[i] = ldm_p0_ptr[i + update_strip_len];
      ldm_p1_ptr[i] = ldm_p1_ptr[i + update_strip_len];
      ldm_v_ptr[i]  = ldm_v_ptr[i+ update_strip_len];
    }

    for (int i = 0; i < update_strip_len; i++) {
      ldm_p0_ptr[max_strip_4_ldm + i] = p0tmp[i];
      ldm_p1_ptr[max_strip_4_ldm + i] = p1tmp[i];
      ldm_v_ptr[max_strip_4_ldm + i]  = vtmp[i];
    }

    for (int i = 0; i < num_strip_to_copy; i++) {
      while (put_reply[i] != 1);
    }
    for (int i = 0; i < next_num_strip_to_copy; i++) {
      while (get_reply_p0[i] != 1);
      while (get_reply_p1[i] != 1);
      while (get_reply_v[i] != 1);
    }


    if (debug) printf("next_num_strip_to_copy: %d\n", next_num_strip_to_copy);
    host_istrip += next_num_strip_to_copy;
    cur_nx = max_strip_4_ldm - update_strip_len + next_num_strip_to_copy;
    num_strip_to_copy = next_num_strip_to_copy;
  } while (num_strip_to_copy > 0);

  /*while (host_istrip < strip_len) {*/
    /*num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);*/

    /*get_reply = 0;*/
    /*for (int i = 0; i < num_strip_to_copy; i++) {*/
      /*athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);*/
      /*athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);*/
      /*athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm-update_strip_len+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);*/
    /*}*/
    /*while (get_reply != 3 * num_strip_to_copy);*/

    /*if (debug) {*/
      /*printf("before stencil\n");*/
      /*print(ldm_p0_ptr, max_strip_4_ldm, nz);*/
    /*}*/

    /*fd4t10s_zjh_2d_vtrans_3(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, ldm_u2_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);*/

    /*if (debug) {*/
      /*printf("after stencil: %d\n", host_istrip);*/
      /*print(ldm_p2_ptr, update_strip_len, nz-2*D);*/
    /*}*/


    /*///copy the data out*/
    /*put_reply = 0;*/
    /*for (int i = 0; i < num_strip_to_copy; i++) {*/
      /*int dst_istrp = dstPos;*/
      /*int src_istrp = i;*/
      /*athread_put(PE_MODE, ldm_p2_ptr[src_istrp], arg.p0+(dst_istrp)*nz + D, (nz-2*D)*sizeof *arg.p0, (void *)&put_reply, 0, 0);*/
      /*dstPos++;*/
    /*}*/
    /*while (put_reply != num_strip_to_copy);*/

    /*/// copy next strip in*/
    /*for (int i = 0; i < update_strip_len; i++) {*/
      /*p0tmp[i] = ldm_p0_ptr[i];*/
      /*p1tmp[i] = ldm_p1_ptr[i];*/
      /*vtmp[i] = ldm_v_ptr[i];*/
    /*}*/
    /*for (int i = 0; i < max_strip_4_ldm - update_strip_len; i++) {*/
      /*ldm_p0_ptr[i] = ldm_p0_ptr[i + update_strip_len];*/
      /*ldm_p1_ptr[i] = ldm_p1_ptr[i + update_strip_len];*/
      /*ldm_v_ptr[i]  = ldm_v_ptr[i+ update_strip_len];*/
    /*}*/

    /*for (int i = 0; i < update_strip_len; i++) {*/
      /*ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i] = p0tmp[i];*/
      /*ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i] = p1tmp[i];*/
      /*ldm_v_ptr[max_strip_4_ldm - update_strip_len + i]  = vtmp[i];*/
    /*}*/

    /*host_istrip += num_strip_to_copy;*/

  /*}*/

  if (debug) {
    printf("tid: %d, exit %s\n", tid, __PRETTY_FUNCTION__);
  }
  if (debug) {
    printf("tid: %d, exit %s\n", tid, __PRETTY_FUNCTION__);
  }

}
