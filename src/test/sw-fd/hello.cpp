#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>
#include <cmath>
#include <cstring>
#include <unistd.h>

#include "hello_arg.h"

extern "C" {
#include "athread.h"
void slave_func(void *);
}

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
    stencil_kernel(p0, p1, v, nx, nz);
    //stencil_kernel_2(p0_ptr, p1_ptr, v_ptr, nx, nz);
//    print(p0_ptr[0], nx, nz);
    int s = std::accumulate(p0, p0 + nx * nz, 0);
    sum.push_back(s);
    std::printf("it: %d, sum of p0: %d\n", it, s);

    std::swap(p0, p1);
  }

  return sum;
}


std::vector<int> run_parallel(int *p0, int *p1, int *v, int nx, int nz, int nt) {
  pthread_t threads[NUM_THREADS];
  struct args_t arg[NUM_THREADS];

  athread_init();

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
     //int rc = pthread_create(&threads[tid], NULL, stencil_wrapper2, (void *)&arg[tid]);
     int rc = __real_athread_create(tid, (void *)slave_func, &arg[tid]);

     if (rc){
       printf("ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
     }
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
    //for (int tid = 0; tid < NUM_THREADS; tid++) {
      //pthread_join(threads[tid], NULL);
    //}
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
  int nx = NX;
  int nz = NZ;
  int nt = NT;
  std::vector<int> p0(nx * nz);
  std::vector<int> p1(nx * nz);
  std::vector<int> v(nx * nz);

  //if(pthread_barrier_init(&barr, NULL, NUM_THREADS))
  //{
    //printf("Could not create a barrier\n");
    //return -1;
  //}

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
