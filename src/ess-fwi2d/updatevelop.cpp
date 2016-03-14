/*
 * updatevelop.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#include "updatevelop.h"

UpdateVelOp::UpdateVelOp() {
  // TODO Auto-generated constructor stub

}

void update_vel(const float *vel, const float *grad, float size, float steplen, float vmin, float vmax, float *new_vel) {
  if (vmax <= vmin) {
    ERROR() << format("vmax(%f) < vmin(%f)") % vmax % vmin;
    exit(0);
  }

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
