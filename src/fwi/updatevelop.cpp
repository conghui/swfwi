/*
 * updatevelop.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#include "updatevelop.h"
#include "logger.h"

static void update_vel(float *new_vel, const float *vel, const float *grad, int size, float steplen, float vmin, float vmax) {
  for (int i = 0; i < size; i++) {
    new_vel[i] = vel[i] + steplen * grad[i];
    if (new_vel[i] > vmax) {
      new_vel[i] = vmax;
    }
    if (new_vel[i] < vmin) {
      new_vel[i] = vmin;
    }
  }
}

UpdateVelOp::UpdateVelOp(float _vmin, float _vmax, float dx, float dt) {
  vmin = (dx / dt / _vmax) * (dx / dt / _vmax);
  vmax = (dx / dt / _vmin) * (dx / dt / _vmin);

  if (vmax <= vmin) {
    ERROR() << format("vmax(%f) < _vmin(%f)") % vmax % vmin;
    exit(0);
  }
}

void UpdateVelOp::update(Velocity& newVel,
    const Velocity& vel, const std::vector<float>& grad,
    float steplen) const {
  update_vel(&newVel.dat[0], &vel.dat[0], &grad[0], newVel.dat.size(), steplen, vmin, vmax);
}
