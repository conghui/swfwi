#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>
#include <cmath>
#include <cstring>

#define NUM_THREADS 3
#define D 2

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

void run_serial(int *p0, int *p1, int *v, int nx, int nz, int nt) {
  for (int it = 0; it < nt; it++) {
    stencil_kernel(p0, p1, v, nx, nz);
    std::printf("it: %d, sum of p0: %d\n", it, std::accumulate(p0, p0 + nx * nz, 0));
    std::swap(p0, p1);
  }
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
    //print(&local_p0[0], (xend - xbegin), arg.nz);

    //printf("\nlocal_p1:\n");
    //print(&local_p1[0], (xend - xbegin), arg.nz);
    /// calculate sum before output

  }
  std::memcpy(arg.p0 + (p0_x_output_begin * arg.nz), &local_p0[local_p0_x_output_begin*arg.nz], local_p0_x_output_len * arg.nz * sizeof(*arg.p0));

  return NULL;
}

void run_parallel(int *p0, int *p1, int *v, int nx, int nz, int nt) {
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

void *PrintHello(void *threadid)
{
   long tid;
   tid = (long)threadid;
   printf("Hello World! It's me, thread #%ld!\n", tid);
   pthread_exit(NULL);
}


int main(int argc, char *argv[])
{
  int nx = 117;
  int nz = 47;
  int nt = 10;
  std::vector<int> p0(nx * nz);
  std::vector<int> p1(nx * nz);
  std::vector<int> v(nx * nz);

  std::generate(p0.begin(), p0.end(), myrand);
  std::generate(p1.begin(), p1.end(), myrand);
  std::fill(v.begin(), v.end(), 2);

  //std::printf("p0:\n");
  //print(&p0[0], nx, nz);

  //std::printf("\np1:\n");
  //print(&p1[0], nx, nz);

  //std::printf("\nv:\n");
  //print(&v[0], nx, nz);

  {
    std::vector<int> tmp_p0 = p0;
    std::vector<int> tmp_p1 = p1;
    std::vector<int> tmp_v = v;
    run_serial(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
  }

  //run_serial(&p0[0], &p1[0], &v[0], nx, nz, nt);
  //std::printf("p0:\n");
  //print(&p0[0], nx, nz);

  {
    std::vector<int> tmp_p0 = p0;
    std::vector<int> tmp_p1 = p1;
    std::vector<int> tmp_v = v;
    run_parallel(&tmp_p0[0], &tmp_p1[0], &tmp_v[0], nx, nz, nt);
  }

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
