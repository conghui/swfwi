/*
 * ipdatesteplenop.h
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_UPDATESTEPLENOP_H_
#define SRC_ESS_FWI2D_UPDATESTEPLENOP_H_

#include <vector>
#include "damp4t10d.h"
#include "updatevelop.h"

class UpdateSteplenOp {
public:
  UpdateSteplenOp(const Damp4t10d &fmMethod, const UpdateVelOp &updateVelOp, int max_iter_select_alpha3, float maxdv);

  void bindEncSrcObs(const std::vector<float> &encsrc, const std::vector<float> &encobs);
  void calsteplen(const std::vector<float> &grad, float obj_val1, int iter, float lambdaX, float lambdaZ, float &steplen, float &objval);

private:
  float calobjval(const std::vector<float> &grad, float steplen) const;
  bool refineAlpha(const std::vector<float> &grad, float obj_val1, float maxAlpha3, float &_alpha2, float &_obj_val2, float &_alpha3, float &_obj_val3) const;
  void initAlpha23(float maxAlpha3, float &initAlpha2, float &initAlpha3);

private:
  struct PreservedAlpha {
    float alpha;
    bool init;

    PreservedAlpha() : alpha(0), init(false) {}
  } preservedAlpha;

private:
  const Damp4t10d &fmMethod;
  const UpdateVelOp &updateVelOp;
  const std::vector<float> *encsrc;
  const std::vector<float> *encobs;

  int max_iter_select_alpha3;
  float maxdv;

  float lambdaX, lambdaZ;
};

#endif /* SRC_ESS_FWI2D_UPDATESTEPLENOP_H_ */
