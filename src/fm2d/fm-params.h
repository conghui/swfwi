/*
 * fm-params.h
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_FM_PARAMS_H_
#define SRC_FM2D_FM_PARAMS_H_


extern "C"
{
#include <rsf.h>
}

class FmParams {
public:
  static FmParams &instance();

private:
  FmParams();
  FmParams(const FmParams &);
  void operator=(const FmParams &);
  ~FmParams();
  void getInputParams();
  void putOutputParams();
  void check();

private:
  static FmParams *ins;

public:
  sf_file vinit;
  sf_file shots;
  int nb;
  int nz;
  int nx;
  float dz;
  float dx;
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
};

#endif /* SRC_FM2D_FM_PARAMS_H_ */
