/*
 * damp4t10d.h
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_DAMP4T10D_H_
#define SRC_FM2D_DAMP4T10D_H_

#include "velocity.h"
#include "shot-position.h"

class Damp4t10d {
public:
  Damp4t10d(float dt, float dx, int nb);
  Velocity expandVelocity(const Velocity &vel);

  void stepForward(float *p0, float *p1);
  void bindVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, const ShotPosition &pos);
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos);

private:
  void fd4t10s_bd1_2d(float *prev_wave, const float *curr_wave);
private:
  const static int FDLEN = 5;

private:
  const Velocity *vel;
  float dt;
  float dx;
  int nb;
};

#endif /* SRC_FM2D_DAMP4T10D_H_ */
