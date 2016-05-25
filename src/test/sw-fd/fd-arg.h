#ifndef ARG_H_NIYPCBT7
#define ARG_H_NIYPCBT7

#define NUM_THREADS 64
#define D 6
#define NX 40960
#define NZ 193 /*193*/
#define NT 1

struct args_t {
  float *p0;
  float *p1;
  float *v;
  int nx;
  int nz;
  int nthreads;
} ;

#endif /* end of include guard: ARG_H_NIYPCBT7 */
