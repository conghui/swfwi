/*
 * velocity.h
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_VELOCITY_H_
#define SRC_FM2D_VELOCITY_H_

#include <vector>

class Velocity {
public:
  Velocity(int _nx, int _nz);
  Velocity (const std::vector<float> &dat, int nx, int nz);
public:
  std::vector<float> dat;
  int nx;
  int nz;
};

#endif /* SRC_FM2D_VELOCITY_H_ */
