/*
 * damp4t10d.h
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_DAMP4T10D_H_
#define SRC_FM2D_DAMP4T10D_H_

#include "velocity.h"
class Damp4t10d {
public:
  Damp4t10d(float dt, float dx);
  Velocity transformVelocityForModeling(const Velocity &vel);

  void stepForward(float *p0, float *p1);
  void setVelocity(const Velocity &_vel);

private:
  const static int FDLEN = 5;

private:
  float dt;
  float dx;
};

#endif /* SRC_FM2D_DAMP4T10D_H_ */
