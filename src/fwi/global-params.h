/*
 * global-params.h
 *
 *  Created on: Feb 23, 2016
 *      Author: rice
 */

#ifndef SRC_FWI_GLOBAL_PARAMS_H_
#define SRC_FWI_GLOBAL_PARAMS_H_

extern "C"
{
#include <rsf.h>
}

class GlobalParams {
public:
  static GlobalParams &instance();

private:
  GlobalParams();
  GlobalParams(const GlobalParams &);
  void operator=(const GlobalParams &);
  ~GlobalParams();
  void getInputParams();
  void putOutputParams();

private:
  static GlobalParams *ins;

public:
  sf_file vinit;        /* initial velocity model, unit=m/s */
  sf_file shots;        /* recorded shots from exact velocity model */
  sf_file vupdates;     /* updated velocity in iterations */
  sf_file grads;        /* gradient in iterations */
  sf_file illums;       /* source illumination in iterations */
  sf_file objs;         /* values of objective function in iterations */
  bool verb;            // verbosity
  int nz;
  int nx;
  float dz;
  float dx;
  bool precon;
  int niter;
  int rbell;
  int nt;
  int ng;
  int ns;
  float dt;
  float amp;
  float fm;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  int jsx;
  int jsz;
  int jgx;
  int jgz;
  int csd;
};

#endif /* SRC_FWI_GLOBAL_PARAMS_H_ */
