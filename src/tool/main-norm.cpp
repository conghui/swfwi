
extern "C" {
#include <rsf.h>
}

#include <cmath>
#include <vector>
#include "environment.h"

namespace {
class Params {
public:
  Params();
  ~Params();

private:
  Params(const Params &);
  void operator=(const Params &);

public:
  sf_file realvel;
  sf_file velset;
  sf_file l1norm;
  sf_file l2norm;
  int nz;
  int nx;
  int nt;
};


Params::Params() {
  /*< set up I/O files >*/
  realvel=sf_input ("realvel");   /* initial velocity model, unit=m/s */
  velset = sf_input("velset");
  l1norm = sf_output("l1norm");
  l2norm = sf_output("l2norm");

  /* get parameters for forward modeling */
  if (!sf_histint(velset,"n1",&nz)) sf_error("no n1 in velset");
  if (!sf_histint(velset,"n2",&nx)) sf_error("no n2 in velset");
  if (!sf_histint(velset,"n3",&nt)) sf_error("no n3 in velset");

  int tnz, tnx;
  if (!sf_histint(realvel, "n1", &tnz)) sf_error("no n1 in realvel");
  if (!sf_histint(realvel, "n2", &tnx)) sf_error("no n2 in realvel");
  if (!(tnz == nz && tnx == nx)) sf_error("n1 n2 of velset and realvel doesnot match");

  sf_putint(l1norm, "n1", nt);
  sf_putint(l1norm, "o1", 1);
  sf_putint(l1norm, "d1", 1);
  sf_putstring(l1norm, "label1", "Iteration");

  sf_putint(l2norm, "o1", 1);
  sf_putint(l2norm, "d1", 1);
  sf_putint(l2norm, "n1", nt);
  sf_putstring(l2norm, "label1", "Iteration");
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

  int nz = params.nz;
  int nx = params.nx;
  int nt = params.nt;
  int velsize = nz * nx;
  std::vector<float> realvel(velsize, 0);
  std::vector<float> vel(velsize, 0);
  std::vector<float> l1norm(nt);
  std::vector<float> l2norm(nt);

  sf_floatread(&realvel[0], realvel.size(), params.realvel);

  for (int it = 0; it < nt; it++) {
    sf_floatread(&vel[0], velsize, params.velset);

    float norm1 = 0;
    float norm2 = 0;
    float sumSquare = 0;
    float sum = 0;
    for (int i = 0; i < velsize; i++) {
      float accvel = 1.0f / realvel[i];     /// slowness
      float curvel = 1.0f / vel[i];      /// slowness
      float absdiff = std::abs(accvel - curvel);

      norm2     += absdiff * absdiff;
      sumSquare += accvel * accvel;

      norm1  += absdiff;
      sum    += accvel;
    }
    l1norm[it] = norm1 / sum;
    l2norm[it] = norm2 / sumSquare;
  }

  sf_floatwrite(&l1norm[0], l1norm.size(), params.l1norm);
  sf_floatwrite(&l2norm[0], l2norm.size(), params.l2norm);

  return 0;
}

