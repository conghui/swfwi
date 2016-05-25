#ifndef HELLO_ARG_H_OGDEAUCT
#define HELLO_ARG_H_OGDEAUCT

#define NUM_THREADS 3
#define D 2
#define NX 19
#define NZ 6
#define NT 1

struct args_t {
  int *p0;
  int *p1;
  int *v;
  int nx;
  int nz;
  int nthreads;
  int tid;
} ;

#endif /* end of include guard: HELLO_ARG_H_OGDEAUCT */
