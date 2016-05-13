#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>
#include <math>

#define NUM_THREADS 3

int myrand() {
  return rand() % 10;
}

/**
 * we use integer instead of floating point numbers to ensure correctness
 * between different methods.
 */
void stencil_kernel(int *p0, const int *p1, const int *v, int nx, int nz) {
  const int D = 2;
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
    std::printf("sum of p0: %d\n", std::accumulate(p0, p0 + nx * nz, 0));
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
  std::printf("nx: %d, nz: %d, tid: %d\n", arg.nx, arg.nz, arg.tid);
  int nstrip_per_thread = std::ceil(1.0 * arg.nx / arg.nthreads);

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
    std::printf("sum of p0: %d\n", std::accumulate(p0, p0 + nx * nz, 0));
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
  int nx = 40;
  int nz = 20;
  int nt = 1;
  std::vector<int> p0(nx * nz);
  std::vector<int> p1(nx * nz);
  std::vector<int> v(nx * nz);

  std::generate(p0.begin(), p0.end(), myrand);
  std::generate(p1.begin(), p1.end(), myrand);
  std::fill(v.begin(), v.end(), 2);

  //run_serial(&p0[0], &p1[0], &v[0], nx, nz, nt);
  run_parallel(&p0[0], &p1[0], &v[0], nx, nz, nt);

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
