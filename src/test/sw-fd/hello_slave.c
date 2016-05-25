#include "slave.h"
#include "hello_arg.h"

#define max_stencil_len  (2 * D + 1)
#define max_strip_4_ldm  (7)
#define update_strip_len  (max_strip_4_ldm - max_stencil_len + 1)

#define MAX(a, b) (a) > (b) ? (a) : (b)
#define MIN(a, b) (a) < (b) ? (a) : (b)

__thread_local int ldm_p0[max_strip_4_ldm + update_strip_len][NZ];
__thread_local int ldm_p1[max_strip_4_ldm + update_strip_len][NZ];
__thread_local int ldm_v[max_strip_4_ldm + update_strip_len][NZ];
__thread_local int ldm_u2[max_strip_4_ldm][NZ];

__thread_local int *ldm_p0_ptr[max_strip_4_ldm + update_strip_len];
__thread_local int *ldm_p1_ptr[max_strip_4_ldm + update_strip_len];
__thread_local int *ldm_v_ptr[max_strip_4_ldm  + update_strip_len];
__thread_local int *ldm_u2_ptr[max_strip_4_ldm];

__thread_local int *p0tmp[update_strip_len];
__thread_local int *p1tmp[update_strip_len];
__thread_local int *vtmp[update_strip_len];

__thread_local int ldm_p2[update_strip_len][NZ-2*D];
__thread_local int *ldm_p2_ptr[update_strip_len];

int sleep(float t) {
  int s = 0;
  for (int i = 0; i < t * 1e9; i++) {
    s++;
  }
  return s;
}

void stencil_kernel_3(int **p0, int **p1, int **p2, int **v, int nx, int nz) {
  int a[D];
  a[0] = 2;
  a[1] = 3;

  for (int ix = D; ix < nx - D; ix++) {
    for (int iz = D; iz < nz - D; iz++) {
      p2[ix][iz] = p0[ix][iz] * v[ix][iz] +
                a[0] * (p1[ix][iz-1] + p1[ix][iz+1] + p1[ix-1][iz] + p1[ix+1][iz]) +
                a[1] * (p1[ix][iz-2] + p1[ix][iz+2] + p1[ix-2][iz] + p1[ix+2][iz]);
    }
  }
}

void stencil_kernel_4(int **p0, int **p1, int **p2, int **v, int nx, int nz) {
  int a[D];
  a[0] = 2;
  a[1] = 3;

  for (int ix = D; ix < nx - D; ix++) {
    for (int iz = D; iz < nz - D; iz++) {
      p2[ix-D][iz-D] = p0[ix][iz] * v[ix][iz] +
                a[0] * (p1[ix][iz-1] + p1[ix][iz+1] + p1[ix-1][iz] + p1[ix+1][iz]) +
                a[1] * (p1[ix][iz-2] + p1[ix][iz+2] + p1[ix-2][iz] + p1[ix+2][iz]);
    }
  }
}


void print2d(int **A, int nx, int nz) {
  for (int iz = 0; iz < nz; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      printf("%4d", A[ix][iz]);
    }
    printf("\n");
  }
}

void print1d(int *A, int nx, int nz) {
  for (int iz = 0; iz < nz; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      printf("%4d", A[ix * nz + iz]);
    }
    printf("\n");
  }
}

void func1(void *_arg) {
  volatile unsigned int get_reply, put_reply;
  struct args_t arg     = *(struct args_t *)_arg;
  int tid               = athread_get_id(-1);
  int nstrip_per_thread = (arg.nx + arg.nthreads - 1) / arg.nthreads; // # of rows for each threads
  int strip_begin       = MAX(nstrip_per_thread * tid - D, 0); // x starts for each threads
  int strip_end         = MIN(nstrip_per_thread * (arg.tid + 1) + D, arg.nx); /// x ends for each threads
  int strip_len         = MAX(strip_end - strip_begin, 0);
  int host_istrip;
  int nz                = NZ;
  /*if (arg.tid == 0) {*/
    //printf("original:\n");

    //for (int ix = 0; ix < max_strip_4_ldm; ix++) {
      //for (int iz = 0; iz < nz; iz++) {
        //printf("%4d", arg.p0[ix * nz + iz]);
      //}
      //printf("\n");
    //}
    //printf("\n");
    //printf("copied:\n");
  /*}*/

  int debug = tid == 0;
  if (debug) printf("strip_begin: %d, strip_end: %d, strip_len: %d\n", strip_begin, strip_end, strip_len);

  if (debug) printf("inside function: %s\n", __PRETTY_FUNCTION__);
  /*if (debug) print1d(&arg.p0[strip_begin*nz], NX, nz);*/

  for (host_istrip = 0; host_istrip < max_strip_4_ldm; host_istrip++) {
    ldm_p0_ptr[host_istrip] = ldm_p0[host_istrip];
    ldm_p1_ptr[host_istrip] = ldm_p1[host_istrip];
    ldm_v_ptr[host_istrip] = ldm_v[host_istrip];
  }

  for (int i = 0; i < update_strip_len; i++) {
    ldm_p2_ptr[i] = ldm_p2[i];
  }

  /*return;*/

  int num_strip_to_copy = MIN(max_strip_4_ldm - update_strip_len, strip_len);
  if (debug) printf("max_strip_4_ldm: %d\n", max_strip_4_ldm);
  if (debug) printf("update_strip_len: %d\n", update_strip_len);
  if (debug) printf("num_strip_to_copy: %d\n", num_strip_to_copy);

  for (host_istrip = 0; host_istrip < num_strip_to_copy; host_istrip++) {
    get_reply = 0;
    athread_get(PE_MODE, arg.p0+(host_istrip+strip_begin)*nz, ldm_p0[host_istrip], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.p1+(host_istrip+strip_begin)*nz, ldm_p1[host_istrip], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.v+(host_istrip+strip_begin)*nz,  ldm_v[host_istrip],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
    while (get_reply != 3);
  }

  /*if(debug) print1d(ldm_p0[0], max_strip_4_ldm, nz);*/
  /*if(debug) print2d(ldm_p0, max_strip_4_ldm, nz);*/
  if(debug) print2d(ldm_p0_ptr, max_strip_4_ldm, nz);

  /*return;*/
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

  int watched_thread = 1;
  while (host_istrip < strip_len) {
    int num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);
    /*if (arg.tid == watched_thread) {*/
      /*printf("before copy\n");*/
      /*print(ldm_p0_ptr, max_strip_4_ldm, nz);*/
    /*}*/

    /// copy the strip
    for (int i = 0; i < num_strip_to_copy; i++) {
      /*std::memcpy(ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i], arg.p0 + (i + host_istrip + strip_begin) * nz, nz * sizeof *arg.p0);*/
      /*std::memcpy(ldm_p2_ptr[max_strip_4_ldm - update_strip_len + i], arg.p0 + (i + host_istrip + strip_begin) * nz, nz * sizeof *arg.p0);*/
      /*std::memcpy(ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i], arg.p1 + (i + host_istrip + strip_begin) * nz, nz * sizeof *arg.p1);*/
      /*std::memcpy(ldm_v_ptr[max_strip_4_ldm - update_strip_len + i], arg.v +   (i + host_istrip + strip_begin ) * nz, nz * sizeof *arg.v);*/
      get_reply = 0;
      athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
      athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
      athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm-update_strip_len+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
      while (get_reply != 3);
    }

    /*if (arg.tid == watched_thread) {*/
      /*printf("before stencil\n");*/
      /*print(ldm_p0_ptr, max_strip_4_ldm, nz);*/
    /*}*/

//    usleep(0.1 * 1e6);
//  pthread_barrier_wait(&barr);
    /*stencil_kernel_3(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);*/
    stencil_kernel_4(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);

    /*if (arg.tid == watched_thread) {*/
      /*printf("after stencil: %d\n", host_istrip);*/
      /*print(ldm_p2_ptr, max_strip_4_ldm, nz);*/
    /*}*/

    //if (host_istrip == 5) {
      //pthread_barrier_wait(&barr);
    //}

    ///copy the data out
    for (int i = 0; i < num_strip_to_copy; i++) {
//      int dst_istrp = host_istrip - D + i;
//      int dst_istrp = host_istrip - update_strip_len + i;
//      int dst_istrp = host_istrip - num_strip_to_copy + i;
      /*int dst_istrp = dstPos;*/
      /*printf("dstpos: %d\n", dstPos);*/
/*//      int src_istrp = max_strip_4_ldm - (D + 1) + i;*/
      /*int src_istrp = D + i;*/
      /*std::memcpy(arg.p0 + (dst_istrp) * nz, ldm_p2_ptr[src_istrp], nz * sizeof *arg.p0);*/

      /*dstPos++;*/

      put_reply = 0;
      int dst_istrp = dstPos;
      /*[>printf("dstpos: %d\n", dstPos);<]*/
      int src_istrp = i;
      /*[>std::memcpy(arg.p0 + (dst_istrp) * nz, ldm_p2_ptr[src_istrp], nz * sizeof *arg.p0);<]*/
      athread_put(PE_MODE, ldm_p2_ptr[src_istrp], arg.p0+(dst_istrp)*nz + D, (nz-2*D) *sizeof *arg.p0, (void *)&put_reply, 0, 0);
      dstPos++;
      while (put_reply != 1);
    }

//    usleep(0.1 * 1e6);
    /*if (arg.tid == watched_thread) {*/
      /*printf("global p:\n");*/
      /*print(arg.p0, nx, nz);*/
    /*}*/

    /// copy next strip in
    /*int *p0tmp[update_strip_len];*/
    /*int *p1tmp[update_strip_len];*/
    /*int *p2tmp[update_strip_len];*/
    /*int *vtmp[update_strip_len];*/
    for (int i = 0; i < update_strip_len; i++) {
      p0tmp[i] = ldm_p0_ptr[i];
      p1tmp[i] = ldm_p1_ptr[i];
      /*p2tmp[i] = ldm_p2_ptr[i];*/
      vtmp[i] = ldm_v_ptr[i];
    }
    for (int i = 0; i < max_strip_4_ldm - update_strip_len; i++) {
      ldm_p0_ptr[i] = ldm_p0_ptr[i + update_strip_len];
      ldm_p1_ptr[i] = ldm_p1_ptr[i + update_strip_len];
      /*ldm_p2_ptr[i] = ldm_p2_ptr[i + update_strip_len];*/
      ldm_v_ptr[i]  = ldm_v_ptr[i+ update_strip_len];
    }

//    for (int i = max_strip_4_ldm - update_strip_len; i < max_strip_4_ldm; i++) {
    for (int i = 0; i < update_strip_len; i++) {
      ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i] = p0tmp[i];
      ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i] = p1tmp[i];
      /*ldm_p2_ptr[max_strip_4_ldm - update_strip_len + i] = p2tmp[i];*/
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

}

void func(void *_arg) {
  volatile unsigned int get_reply, put_reply;
  struct args_t arg     = *(struct args_t *)_arg;
  int tid               = athread_get_id(-1);
  int nstrip_per_thread = (arg.nx + arg.nthreads - 1) / arg.nthreads; // # of rows for each threads
  int strip_begin       = MAX(nstrip_per_thread * tid - D, 0); // x starts for each threads
  int strip_end         = MIN(nstrip_per_thread * (arg.tid + 1) + D, arg.nx); /// x ends for each threads
  int strip_len         = MAX(strip_end - strip_begin, 0);
  int host_istrip;
  int nz                = NZ;

  int debug = tid == 0;
  if (debug) printf("strip_begin: %d, strip_end: %d, strip_len: %d\n", strip_begin, strip_end, strip_len);

  if (debug) printf("inside function: %s\n", __PRETTY_FUNCTION__);

  /* set up the pointer */
  for (int i= 0; i < max_strip_4_ldm + update_strip_len; i++) {
    ldm_p0_ptr[i] = ldm_p0[i];
    ldm_p1_ptr[i] = ldm_p1[i];
    ldm_v_ptr[i] = ldm_v[i];
  }

  for (int i = 0; i < update_strip_len; i++) {
    ldm_p2_ptr[i] = ldm_p2[i];
  }

  /**
   * First, we should prepare enough data for stencil
   */
  int num_strip_to_copy = MIN(max_strip_4_ldm - update_strip_len, strip_len);
  if (debug) printf("max_strip_4_ldm: %d\n", max_strip_4_ldm);
  if (debug) printf("update_strip_len: %d\n", update_strip_len);
  if (debug) printf("num_strip_to_copy: %d\n", num_strip_to_copy);

  for (host_istrip = 0; host_istrip < num_strip_to_copy; host_istrip++) {
    get_reply = 0;
    athread_get(PE_MODE, arg.p0+(host_istrip+strip_begin)*nz, ldm_p0[host_istrip], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.p1+(host_istrip+strip_begin)*nz, ldm_p1[host_istrip], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.v+(host_istrip+strip_begin)*nz,  ldm_v[host_istrip],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
    while (get_reply != 3);
  }


  num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);
  for (int i = 0; i < num_strip_to_copy; i++) {
    get_reply = 0;
    athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm-update_strip_len+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
    athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm-update_strip_len+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
    while (get_reply != 3);
  }
  if (debug) printf("num_strip_to_copy: %d\n", num_strip_to_copy);
  if(debug) print2d(ldm_p0_ptr, max_strip_4_ldm, nz);


  host_istrip += num_strip_to_copy;
  int dstPos = strip_begin + D;
  int cur_nx = host_istrip;

  do {
    int next_num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);
    if (debug) printf("host_istrip: %d\n", host_istrip);
    get_reply = 0;
    for (int i = 0; i < next_num_strip_to_copy; i++) {
      athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);
      athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);
      athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);
    }

    /*if (debug) {*/
      /*printf("before stencil:\n");*/
      /*print2d(ldm_p0_ptr, cur_nx, nz);*/
    /*}*/

    /*if (debug) printf("cur_nx: %d\n", cur_nx);*/
    stencil_kernel_4(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, cur_nx, nz);

    ///copy the data out
    for (int i = 0; i < num_strip_to_copy; i++) {
      put_reply = 0;
      int dst_istrp = dstPos;
      int src_istrp = i;
      athread_put(PE_MODE, ldm_p2_ptr[src_istrp], arg.p0+(dst_istrp)*nz + D, (nz-2*D) *sizeof *arg.p0, (void *)&put_reply, 0, 0);
      dstPos++;
      while (put_reply != 1);
    }

    /*if (debug) {*/
      /*printf("global p0:\n");*/
      /*print1d(arg.p0, NX, nz);*/
    /*}*/

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

    while (get_reply != 3 * next_num_strip_to_copy);

    if (debug) printf("next_num_strip_to_copy: %d\n", next_num_strip_to_copy);
    host_istrip += next_num_strip_to_copy;
    cur_nx = max_strip_4_ldm - update_strip_len + next_num_strip_to_copy;
    num_strip_to_copy = next_num_strip_to_copy;

  } while (num_strip_to_copy > 0);


  /*while (host_istrip < strip_len) {*/
    /*num_strip_to_copy = MIN(update_strip_len, strip_len - host_istrip);*/

    /*/// prepare data for next iteration*/
    /*get_reply = 0;*/
    /*for (int i = 0; i < num_strip_to_copy; i++) {*/
      /*athread_get(PE_MODE, arg.p0+(i+host_istrip+strip_begin)*nz, ldm_p0_ptr[max_strip_4_ldm+i], nz*sizeof *arg.p0, (void *)&get_reply, 0, 0, 0);*/
      /*athread_get(PE_MODE, arg.p1+(i+host_istrip+strip_begin)*nz, ldm_p1_ptr[max_strip_4_ldm+i], nz*sizeof *arg.p1, (void *)&get_reply, 0, 0, 0);*/
      /*athread_get(PE_MODE, arg.v+(i+host_istrip+strip_begin)*nz,  ldm_v_ptr[max_strip_4_ldm+i],  nz*sizeof *arg.v,  (void *)&get_reply, 0, 0, 0);*/
    /*}*/

    /*stencil_kernel_4(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);*/

    /*///copy the data out*/
    /*for (int i = 0; i < num_strip_to_copy; i++) {*/
      /*put_reply = 0;*/
      /*int dst_istrp = dstPos;*/
      /*int src_istrp = i;*/
      /*athread_put(PE_MODE, ldm_p2_ptr[src_istrp], arg.p0+(dst_istrp)*nz + D, (nz-2*D) *sizeof *arg.p0, (void *)&put_reply, 0, 0);*/
      /*dstPos++;*/
      /*while (put_reply != 1);*/
    /*}*/

    /*while (get_reply != 3 * num_strip_to_copy);*/

    /*for (int i = 0; i < update_strip_len; i++) {*/
      /*p0tmp[i] = ldm_p0_ptr[i];*/
      /*p1tmp[i] = ldm_p1_ptr[i];*/
      /*[>p2tmp[i] = ldm_p2_ptr[i];<]*/
      /*vtmp[i] = ldm_v_ptr[i];*/
    /*}*/
    /*for (int i = 0; i < max_strip_4_ldm - update_strip_len; i++) {*/
      /*ldm_p0_ptr[i] = ldm_p0_ptr[i + update_strip_len];*/
      /*ldm_p1_ptr[i] = ldm_p1_ptr[i + update_strip_len];*/
      /*[>ldm_p2_ptr[i] = ldm_p2_ptr[i + update_strip_len];<]*/
      /*ldm_v_ptr[i]  = ldm_v_ptr[i+ update_strip_len];*/
    /*}*/

    /*for (int i = 0; i < update_strip_len; i++) {*/
      /*ldm_p0_ptr[max_strip_4_ldm - update_strip_len + i] = p0tmp[i];*/
      /*ldm_p1_ptr[max_strip_4_ldm - update_strip_len + i] = p1tmp[i];*/
      /*[>ldm_p2_ptr[max_strip_4_ldm - update_strip_len + i] = p2tmp[i];<]*/
      /*ldm_v_ptr[max_strip_4_ldm - update_strip_len + i]  = vtmp[i];*/
    /*}*/

    /*host_istrip += num_strip_to_copy;*/

  /*}*/

}
