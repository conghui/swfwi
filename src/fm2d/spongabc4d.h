/*
 * spongabc4d.h
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_SPONGABC4D_H_
#define SRC_FM2D_SPONGABC4D_H_

#include "velocity.h"
#include "fm-params.h"

class SpongAbc4d {
public:
  SpongAbc4d(const Velocity &_vel, const FmParams &_params);
  void stepForward(float *p0, float *p1);

private:
  void applySponge(float *p);
  void transVel();
  void initCoeff();
  void initbndr();

public:
  Velocity vel;
  const FmParams &params;
  std::vector<float> bndr;
  float c0, c11, c12, c21, c22;
};

#endif /* SRC_FM2D_SPONGABC4D_H_ */
