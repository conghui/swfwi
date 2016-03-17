/*
 * i-modeling.cpp
 *
 *  Created on: Mar 2, 2016
 *      Author: rice
 */

#include <cstdlib>
#include "i-modeling.h"
#include "logger.h"

IModeling::IModeling() : vel(NULL)
{

}

IModeling::~IModeling() {
}

void IModeling::bindVelocity(const Velocity& vel) {
  this->vel = &vel;
}
