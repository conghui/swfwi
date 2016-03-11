/*
 * sfutil.cpp
 *
 *  Created on: Mar 3, 2016
 *      Author: rice
 */


extern "C"{
#include <rsf.h>
}

#include <vector>
#include "sfutil.h"

void sfFloatWrite1d(const char *fn, const float *dat, int n1, float d1, float o1) {
  sf_file f = sf_output(fn);
  sf_putint(f, "n1", n1);
  sf_putfloat(f, "d1", d1);
  sf_putfloat(f, "o1", o1);
  sf_floatwrite(const_cast<float *>(dat), n1, f);
}

void sfFloatWrite2d(const char *fn, const float *dat, int n1, int n2,
    float d1, float d2, float o1, float o2) {
  sf_file f = sf_output(fn);
  sf_putint(f, "n1", n1);
  sf_putint(f, "n2", n2);
  sf_putfloat(f, "d1", d1);
  sf_putfloat(f, "d2", d2);
  sf_putfloat(f, "o1", o1);
  sf_putfloat(f, "o2", o2);
  sf_floatwrite(const_cast<float *>(dat), n1 * n2, f);
}

void sfDoubleWrite2d(const char *fn, const double *dat, int n1, int n2,
    float d1, float d2, float o1, float o2) {
  sf_file f = sf_output(fn);
  sf_putint(f, "n1", n1);
  sf_putint(f, "n2", n2);
  sf_putfloat(f, "d1", d1);
  sf_putfloat(f, "d2", d2);
  sf_putfloat(f, "o1", o1);
  sf_putfloat(f, "o2", o2);

  std::vector<float> tmp(dat, dat + n1 * n2);
  sf_floatwrite(&tmp[0], n1 * n2, f);
}
