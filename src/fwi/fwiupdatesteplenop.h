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
#include "fwiupdatevelop.h"

class FwiUpdateSteplenOp {
public:
  FwiUpdateSteplenOp(const Damp4t10d &fmMethod, const FwiUpdateVelOp &updateVelOp, int max_iter_select_alpha3, float maxdv);

  void bindEncSrcObs(const std::vector<float> &encsrc, const std::vector<float> &encobs);
  void calsteplen(const std::vector<float> &grad, float obj_val1, int iter, float &steplen, float &objval);
	void parabola_fit(float alpha1, float alpha2, float alpha3, float obj_val1, float obj_val2, float obj_val3, float maxAlpha3, bool toParabolic, int iter, float &steplen, float &objval);

public:
	float alpha1, alpha2, alpha3, obj_val1, obj_val2, obj_val3;
	float maxAlpha3;

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
  const FwiUpdateVelOp &updateVelOp;
  const std::vector<float> *encsrc;
  const std::vector<float> *encobs;

  int max_iter_select_alpha3;
  float maxdv;
};

#endif /* SRC_ESS_FWI2D_UPDATESTEPLENOP_H_ */
