/*
 * velocity.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include "velocity.h"
#include "logger.h"

Velocity::Velocity(int _nx, int _nz) : dat(_nx *_nz), nx(_nx), nz(_nz) {
}
