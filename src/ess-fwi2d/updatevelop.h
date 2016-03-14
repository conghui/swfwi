/*
 * updatevelop.h
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_UPDATEVELOP_H_
#define SRC_ESS_FWI2D_UPDATEVELOP_H_

class UpdateVelOp {
public:
  UpdateVelOp(float vmin, float vmax);

private:
  float vmin;
  float vmax;
};

#endif /* SRC_ESS_FWI2D_UPDATEVELOP_H_ */
