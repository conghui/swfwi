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
#include "shot-position.h"

class SpongAbc4d {
public:
  SpongAbc4d(float _dt, float _dx, float _dz, int _nb);
  Velocity transformVelocityForModeling(const Velocity &v0);
  void stepForward(float *p0, float *p1);
  void setVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, int ns, const int *sxz, int snz);
  void addSource(float *p, const float *source, const ShotPosition &pos);

private:
  void applySponge(float *p);
  void initCoeff();
  void initbndr();

private:
  const Velocity *vel;
  std::vector<float> bndr;
  float dt, dx, dz;
  int nb;
  float c0, c11, c12, c21, c22;
};

#endif /* SRC_FM2D_SPONGABC4D_H_ */
