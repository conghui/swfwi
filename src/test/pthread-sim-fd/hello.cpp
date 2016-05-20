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
#define D 2
const int nt = 2;
const int nx = 11;
const int nz = 6;
const int max_stencil_len = 2 * D + 1;
const int max_strip_4_ldm = 7; /// must be >= max_stencil_len
const int update_strip_len = max_strip_4_ldm - max_stencil_len + 1;

pthread_barrier_t barr;

int myrand() {
  return rand() % 10;
}

template <typename T>
void print(const T *A, int nx, int nz) {
  for (int iz = 0; iz < nz ; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      std::printf("%4d", A[ix * nz + iz]);
    }
    std::printf("\n");
  }
    std::printf("\n");
}

template <typename T>
void print(T **A, int nx, int nz) {
  for (int iz = 0; iz < nz ; iz++) {
    for (int ix = 0; ix < nx; ix++) {
      std::printf("%4d", A[ix][iz]);
    }
    std::printf("\n");
  }
    std::printf("\n");
}

/**
 * we use integer instead of floating point numbers to ensure correctness
 * between different methods.
 */
void stencil_kernel(int *p0, const int *p1, const int *v, int nx, int nz) {
  int a[D];
  a[0] = 2;
  a[1] = 3;

  for (int ix = D; ix < nx - D; ix++) {
    for (int iz = D; iz < nz - D; iz++) {
      int idx = ix * nz + iz;
      p0[idx] = p0[idx] * v[idx] +
                a[0] * (p1[idx - 1] + p1[idx + 1] + p1[idx - 1 * nz] + p1[idx + 1 * nz]) +
                a[1] * (p1[idx - 2] + p1[idx + 2] + p1[idx - 2 * nz] + p1[idx + 2 * nz]);
    }
  }
}

void stencil_kernel_2(int **p0, int **p1, int **v, int nx, int nz) {
  int a[D];
  a[0] = 2;
  a[1] = 3;

  for (int ix = D; ix < nx - D; ix++) {
    for (int iz = D; iz < nz - D; iz++) {
      p0[ix][iz] = p0[ix][iz] * v[ix][iz] +
                a[0] * (p1[ix][iz-1] + p1[ix][iz+1] + p1[ix-1][iz] + p1[ix+1][iz]) +
                a[1] * (p1[ix][iz-2] + p1[ix][iz+2] + p1[ix-2][iz] + p1[ix+2][iz]);
    }
  }
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

std::vector<int> run_serial(int *p0, int *p1, int *v, int nx, int nz, int nt) {
  int *p0_ptr[nx];
  int *p1_ptr[nx];
  int *v_ptr[nx];

  std::vector<int> sum;
  for (int it = 0; it < nt; it++) {
    for (int ix = 0; ix < nx; ix++) {
      p0_ptr[ix] = &p0[ix * nz];
      p1_ptr[ix] = &p1[ix * nz];
      v_ptr[ix] = &v[ix * nz];
    }
    //stencil_kernel(p0, p1, v, nx, nz);
    stencil_kernel_2(p0_ptr, p1_ptr, v_ptr, nx, nz);
//    print(p0_ptr[0], nx, nz);
    int s = std::accumulate(p0, p0 + nx * nz, 0);
    sum.push_back(s);
    std::printf("it: %d, sum of p0: %d\n", it, s);

    std::swap(p0, p1);
  }

  return sum;
}

struct args_t {
  int *p0;
  int *p1;
  int *v;
  int nx;
  int nz;
  int nthreads;
  int tid;
};

void *stencil_wrapper(void *_arg) {
  struct args_t arg = *(struct args_t *)_arg;
  int nstrip_per_thread = std::ceil(1.0 * arg.nx / arg.nthreads); // # of rows for each threads
  int strip_begin = std::max(nstrip_per_thread * arg.tid - D, 0); // x starts for each threads
  int strip_end =  std::min(nstrip_per_thread * (arg.tid + 1) + D, arg.nx); /// x ends for each threads
  int strip_len = strip_end - strip_begin;

  std::vector<int> local_p0(arg.nz * strip_len);
  std::vector<int> local_p1(arg.nz * strip_len);
  std::vector<int> local_v(arg.nz * strip_len);

  std::memcpy(&local_p0[0], arg.p0 + strip_begin * arg.nz, strip_len * arg.nz * sizeof *arg.p0);
  std::memcpy(&local_p1[0], arg.p1 + strip_begin * arg.nz, strip_len * arg.nz * sizeof *arg.p1);
  std::memcpy(&local_v[0], arg.v + strip_begin * arg.nz, strip_len * arg.nz * sizeof *arg.v);

  if (arg.tid == 0) {
    //printf("\nlocal_p0\n");
    //print(&local_p0[0], (strip_end - strip_begin), arg.nz);
  }
  pthread_barrier_wait(&barr);

  stencil_kernel(&local_p0[0], &local_p1[0], &local_v[0], strip_end - strip_begin, arg.nz);

  int p0_x_output_begin = arg.tid == 0 ? D : nstrip_per_thread * arg.tid;
  int local_p0_x_output_begin = D;
  int local_p0_x_output_len = strip_end - strip_begin - 2 * D;

  if (arg.tid == 0) {
    //std::printf("nstrip_per_thread: %d, xbegin: %d, xend: %d\n", nstrip_per_thread, xbegin, xend);
    //std::printf("tid: %d, p0_x_output_begin: %d, local_p0_x_output_begin: %d, local_p0_x_output_len: %d \n",
        //arg.tid, p0_x_output_begin, local_p0_x_output_begin, local_p0_x_output_len);

    //printf("\nglobal_p0:\n");
    //print(arg.p0, arg.nx, arg.nz);

    //printf("\nlocal_p0\n");
    //print(&local_p0[0], (strip_end - strip_begin), arg.nz);

    //printf("\nlocal_p1:\n");
    //print(&local_p1[0], (xend - xbegin), arg.nz);
    /// calculate sum before output

  }
  std::memcpy(arg.p0 + (p0_x_output_begin * arg.nz), &local_p0[local_p0_x_output_begin*arg.nz], local_p0_x_output_len * arg.nz * sizeof(*arg.p0));

  return NULL;
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

  int ldm_p0[max_strip_4_ldm][nz] = {0};
  int ldm_p1[max_strip_4_ldm][nz] = {0};
  int ldm_p2[max_strip_4_ldm][nz] = {0};
  int ldm_v[max_strip_4_ldm][nz]  = {0};

  int *ldm_p0_ptr[max_strip_4_ldm];
  int *ldm_p1_ptr[max_strip_4_ldm];
  int *ldm_p2_ptr[max_strip_4_ldm];
  int *ldm_v_ptr[max_strip_4_ldm];

  for (host_istrip = 0; host_istrip < max_strip_4_ldm; host_istrip++) {
    ldm_p0_ptr[host_istrip] = ldm_p0[host_istrip];
    ldm_p1_ptr[host_istrip] = ldm_p1[host_istrip];
    ldm_p2_ptr[host_istrip] = ldm_p2[host_istrip];
    ldm_v_ptr[host_istrip] = ldm_v[host_istrip];
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

  int watched_thread = 1;
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
    stencil_kernel_3(ldm_p0_ptr, ldm_p1_ptr, ldm_p2_ptr, ldm_v_ptr, max_strip_4_ldm - (update_strip_len - num_strip_to_copy), nz);

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
    int *p0tmp[update_strip_len];
    int *p1tmp[update_strip_len];
    int *p2tmp[update_strip_len];
    int *vtmp[update_strip_len];
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

void run_parallel0(int *p0, int *p1, int *v, int nx, int nz, int nt) {
  pthread_t threads[NUM_THREADS];
  struct args_t arg[NUM_THREADS];

  for (int it = 0; it < nt; it++) {

    for (int tid = 0; tid < NUM_THREADS; tid++) {
      arg[tid].p0 = p0;
      arg[tid].p1 = p1;
      arg[tid].v = v;
      arg[tid].nx = nx;
      arg[tid].nz = nz;
      arg[tid].nthreads = NUM_THREADS;
      arg[tid].tid = tid;
     int rc = pthread_create(&threads[tid], NULL, stencil_wrapper, (void *)&arg[tid]);

     if (rc){
       printf("ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
     }
    }

    for (int tid = 0; tid < NUM_THREADS; tid++) {
      pthread_join(threads[tid], NULL);
    }
    std::printf("it: %d, sum of p0: %d\n", it, std::accumulate(p0, p0 + nx * nz, 0));
    std::swap(p0, p1);
  }
}

std::vector<int> run_parallel(int *p0, int *p1, int *v, int nx, int nz, int nt) {
  pthread_t threads[NUM_THREADS];
  struct args_t arg[NUM_THREADS];

  std::vector<int> sum;
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
    int s = std::accumulate(p0, p0 + nx * nz, 0);
    sum.push_back(s);
    std::printf("it: %d, sum of p0: %d\n", it, s);
    std::swap(p0, p1);
  }

  return sum;
}

void *PrintHello(void *threadid)
{
   long tid;
   tid = (long)threadid;
   printf("Hello World! It's me, thread #%ld!\n", tid);
   pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
  std::vector<int> p0(nx * nz);
  std::vector<int> p1(nx * nz);
  std::vector<int> v(nx * nz);

  if(pthread_barrier_init(&barr, NULL, NUM_THREADS))
  {
    printf("Could not create a barrier\n");
    return -1;
  }

  std::generate(p0.begin(), p0.end(), myrand);
  std::generate(p1.begin(), p1.end(), myrand);
  std::fill(v.begin(), v.end(), 2);

  std::printf("p0:\n");
  print(&p0[0], nx, nz);

  //std::printf("\np1:\n");
  //print(&p1[0], nx, nz);

  //std::printf("\nv:\n");
  //print(&v[0], nx, nz);

  std::vector<int> sum_serial;
  {
    std::vector<int> tmp_p0 = p0;
    std::vector<int> tmp_p1 = p1;
    std::vector<int> tmp_v = v;
    sum_serial = run_serial(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
    print(&tmp_p0[0], nx, nz);
  }

  //run_serial(&p0[0], &p1[0], &v[0], nx, nz, nt);
  //std::printf("p0:\n");
  //print(&p0[0], nx, nz);

//  {
//    std::vector<int> tmp_p0 = p0;
//    std::vector<int> tmp_p1 = p1;
//    std::vector<int> tmp_v = v;
//    run_parallel0(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
//  }

  std::vector<int> sum_parallel;
  {
    std::vector<int> tmp_p0 = p0;
    std::vector<int> tmp_p1 = p1;
    std::vector<int> tmp_v = v;
    sum_parallel = run_parallel(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
//    print(&tmp_p0[0], nx, nz);
  }

  for (int i = 0; i < nt; i++) {
    if (sum_serial[i] != sum_parallel[i]) {
      printf("error\n");
      exit(0);
    }
  }

  printf("test pass!\n");

  return 0;
   //pthread_t threads[NUM_THREADS];
   //int rc;
   //long t;
   //for(t=0;t<NUM_THREADS;t++){
     //printf("In main: creating thread %ld\n", t);
     //rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);
     //if (rc){
       //printf("ERROR; return code from pthread_create() is %d\n", rc);
       //exit(-1);
       //}
     //}

   //[> Last thing that main() should do <]
   //pthread_exit(NULL);
}
