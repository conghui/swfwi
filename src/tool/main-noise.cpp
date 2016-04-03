
extern "C" {
#include <rsf.h>
}

#include <boost/bind.hpp>
#include <cmath>
#include <cstdio>
#include "common.h"
#include "sum.h"
#include "environment.h"

namespace {
float rand01f() {
  return (float)std::rand() / RAND_MAX;
}

float addNoise(float val, float maxabs, float factor) {
  return val + maxabs * rand01f() * factor;
}

class Params {
public:
  Params();
  ~Params();

private:
  Params(const Params &);
  void operator=(const Params &);

public:
  sf_file shots;
  sf_file shotsnoise;
  int nt;
  int ng;
  int ns;
  float fac;
};


Params::Params() {
  /*< set up I/O files >*/
  shots = sf_input("in");
  shotsnoise = sf_output("out");
  if (!sf_getfloat("fac", &fac)) sf_error("no fac");
  if (!sf_histint(shots,"n1",&nt)) sf_error("no n1 in velset");
  if (!sf_histint(shots,"n2",&ng)) sf_error("no n2 in velset");
  if (!sf_histint(shots,"n3",&ns)) sf_error("no n3 in velset");
}

Params::~Params() {
  sf_close();
}

} /// end of name space

int main(int argc, char* argv[]) {
  /* initialize Madagascar */
  sf_init(argc,argv);
  Environment::setDatapath();

  Params params;

  int nt = params.nt;
  int ng = params.ng;
  int ns = params.ns;

  std::vector<float> shot(ng * nt, 0);
  for (int ishot = 0; ishot < ns; ishot++) {
    sf_floatread(&shot[0], shot.size(), params.shots);
    float maxabs = *std::max_element(shot.begin(), shot.end(), abs_less<float>);
    std::transform(shot.begin(), shot.end(), shot.begin(), boost::bind(addNoise, _1, maxabs, params.fac));
    sf_floatwrite(&shot[0], shot.size(), params.shotsnoise);
  }

  return 0;
}

