/*
 * i-modeling.h
 *
 *  Created on: Mar 2, 2016
 *      Author: rice
 */

#ifndef SRC_MODELING_I_MODELING_H_
#define SRC_MODELING_I_MODELING_H_

#include "velocity.h"

class IModeling {
public:
  IModeling();
  void bindVelocity(const Velocity &vel);

public: // virtual
  virtual Velocity expandDomain(const Velocity &v0) const = 0;
  virtual ~IModeling();

protected:
  const Velocity *vel;
};

#endif /* SRC_MODELING_I_MODELING_H_ */
