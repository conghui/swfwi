#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <numeric>
#include <functional>
#include <sys/time.h>

extern "C" {
#include "add.h"
}

float myrand() {
  return std::rand() * 1.0 / RAND_MAX;
}

void test()
{
  int nx = 533;
  int nz = 193;
  int nb = 36;

  std::vector<float> prev(nx * nz);
  std::vector<float> curr(nx * nz);
  std::vector<float> vel(nx * nz);
  std::vector<float> u2(nx * nz);

  std::generate(curr.begin(), curr.end(), myrand);
  std::fill(vel.begin(), vel.end(), 1500);

  struct timeval start;
  struct timeval stop;

  gettimeofday(&start, NULL);
  fd4t10s_damp_zjh_2d_vtrans(&prev[0], &curr[0], &vel[0], &u2[0], nx, nz, nb);
  gettimeofday(&stop, NULL);

  double elapsed = (stop.tv_sec - start.tv_sec) * 10e6 + (stop.tv_usec - start.tv_usec);
  elapsed /= 1000;

  std::cout << "elapsed: " << elapsed << "ms\n";
  std::cout << "sum prev: " << std::accumulate(prev.begin(), prev.end(), 0.0f) << std::endl;

}
