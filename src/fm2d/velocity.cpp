/*
 * velocity.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include "velocity.h"
#include "logger.h"

Velocity::Velocity(const std::vector<float> &vel, int _nx, int _nz) :
  dat(vel), nx(_nx), nz(_nz) {
  if (dat.size() != size_t(nx * nz)) {
    ERROR() << format("%d elements in velocity, but nx*nz=%d") % dat.size() % (nx * nz);
    exit(0);
  }
}

Velocity::Velocity(int _nx, int _nz) : dat(_nx *_nz), nx(_nx), nz(_nz) {
}
