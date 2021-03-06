/*
 * essfwiframework.h
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_
#define SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_

#include "damp4t10d.h"
#include "fwiupdatevelop.h"
#include "fwiupdatesteplenop.h"
#include "random-code.h"

class FwiFramework {
public:
  FwiFramework(Damp4t10d &fmMethod, const FwiUpdateSteplenOp &updateSteplenOp,
                  const FwiUpdateVelOp &updateVelOp, const std::vector<float> &wlt,
                  const std::vector<float> &dobs);

  void epoch(int iter);
  void writeVel(sf_file file) const;
  float getUpdateObj() const;
  float getInitObj() const;

private:
  static const int ESS_SEED = 1;

private:
  Damp4t10d &fmMethod;
  FwiUpdateSteplenOp updateStenlelOp;
  const FwiUpdateVelOp &updateVelOp;
  const std::vector<float> &wlt;  /// wavelet
  const std::vector<float> &dobs; /// actual observed data (nt*ng*ns)
  RandomCodes essRandomCodes;

private: /// propagate from other construction
  int ns;
  int ng;
  int nt;
  int nx;
  int nz;
  float dx;
  float dt;

private:
  std::vector<float> g0;               /// gradient in previous step
  std::vector<float> updateDirection;
  float updateobj;
  float initobj;
	float obj_val4;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
