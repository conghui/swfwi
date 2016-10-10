/*
 * ipdatesteplenop.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#include <set>
#include <cmath>
#include "updatesteplenop.h"
#include "logger.h"
#include "common.h"
#include "parabola-vertex.h"
#include "sum.h"
#include "ReguFactor.h"

namespace {
typedef std::pair<float, float> ParaPoint;

bool parabolicLessComp(const ParaPoint &a, const ParaPoint &b) {
  return a.second - b.second < 1e-10;
}



void calMaxAlpha2_3(const Velocity &exvel,  const float *grad, float dt, float dx, float maxdv,
                    float &ret_alpha2, float &ret_alpha3) {
  const int nx = exvel.nx;
  const int nz = exvel.nz;

  const std::vector<float> &vel = exvel.dat;
  float alpha2 = FLT_MAX;
  for (int i = 0; i < nx * nz; i++) {
    float tmpv = dx / (dt * std::sqrt(vel[i]));
    tmpv -= maxdv;
    tmpv = (dx / (dt * tmpv)) * (dx / (dt * tmpv));
    if (std::fabs(grad[i]) < 1e-10 ) {
      continue;
    }
    if (alpha2 > (tmpv - vel[i]) / std::fabs(grad[i])) {
      alpha2 = (tmpv - vel[i]) / std::fabs(grad[i]);
    }
  }

  /// return the value
  ret_alpha2 = alpha2;
  ret_alpha3 =  2 * alpha2;
}

} /// end of name space

UpdateSteplenOp::UpdateSteplenOp(const Damp4t10d &fmMethod, const UpdateVelOp &updateVelOp,
    int max_iter_select_alpha3, float maxdv) :
  fmMethod(fmMethod), updateVelOp(updateVelOp), encsrc(NULL), encobs(NULL),
  max_iter_select_alpha3(max_iter_select_alpha3), maxdv(maxdv)
{

}

float UpdateSteplenOp::calobjval(const std::vector<float>& grad,
    float steplen) const {
  int nx = fmMethod.getnx();
  int nz = fmMethod.getnz();
  int nt = fmMethod.getnt();

  const Velocity &oldVel = fmMethod.getVelocity();
  Velocity newVel(nx, nz);
  updateVelOp.update(newVel, oldVel, grad, steplen);

  Damp4t10d updateMethod = fmMethod;
  updateMethod.bindVelocity(newVel);

  //forward modeling
  int ng = fmMethod.getng();
  std::vector<float> dcal(nt * ng);
  updateMethod.EssForwardModeling(*encsrc, dcal);

  updateMethod.removeDirectArrival(&dcal[0]);

  std::vector<float> vdiff(nt * ng, 0);
  vectorMinus(*encobs, dcal, vdiff);
  float val = cal_objective(&vdiff[0], vdiff.size());

    if (!(lambdaX == 0 && lambdaZ == 0)) {
      ReguFactor fac(&newVel.dat[0], nx, nz, lambdaX, lambdaZ);
      val += fac.getReguTerm();
    }

  DEBUG() << format("curr_alpha = %e, pure object value = %e") % steplen % val;

  return val;
}

bool UpdateSteplenOp::refineAlpha(const std::vector<float> &grad, float obj_val1, float maxAlpha3,
    float& _alpha2, float& _obj_val2, float& _alpha3, float& _obj_val3) const {

  TRACE() << "SELECTING THE RIGHT OBJECTIVE VALUE 3";

  float alpha3 = _alpha3;
  float alpha2 = _alpha2;
  float obj_val2, obj_val3 = 0;

  obj_val2 = calobjval(grad, alpha2);
  obj_val3 = calobjval(grad, alpha3);

  DEBUG() << "BEFORE TUNNING";
  DEBUG() << __FUNCTION__ << format(" alpha1 = %e, obj_val1 = %e") % 0. % obj_val1;
  DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % alpha2 % obj_val2;
  DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % alpha3 % obj_val3;


  TRACE() << "maintain a set to store alpha2 that we ever tuned";
  std::set<ParaPoint, bool (*)(const ParaPoint &, const ParaPoint &) > tunedAlpha(parabolicLessComp);
  tunedAlpha.insert(std::make_pair(alpha2, obj_val2));

  DEBUG() << "BEGIN TUNING";
  /// obj_val2 might be quite large, so we should make it smaller by halfing alpha2

  int iter = 0;
  for (; iter < max_iter_select_alpha3 && obj_val2 > obj_val1; iter++) {

    /// pass the property of alpha2 to alpha3
    alpha3 = alpha2;
    obj_val3 = obj_val2;

    /// update alpha2
    alpha2 /= 2;
    obj_val2 = calobjval(grad, alpha2);

    /// store it
    tunedAlpha.insert(std::make_pair(alpha2, obj_val2));
    DEBUG() << __FUNCTION__ << format(" iter = %d, alpha2 = %e, obj_val2 = %e") % iter % alpha2 % obj_val2;
    DEBUG() << __FUNCTION__ << format(" iter = %d, alpha3 = %e, obj_val3 = %e\n") % iter % alpha3 % obj_val3;
  }

  DEBUG() << "SELECT A BETTER ALPHA2 IN " << iter << " ITERS";


  DEBUG() << "tunedAlpha size: " << tunedAlpha.size();
  for (std::set<ParaPoint, bool (*)(const ParaPoint &, const ParaPoint &) >::iterator it = tunedAlpha.begin();
      it != tunedAlpha.end(); ++it) {
    DEBUG() << format("alpha %e, obj %e") % it->first % it->second;
  }

  TRACE() << "check if we need to forward tuning";
  TRACE() << "after halfing in the previous step, obj_val2 might still be larger than obj_val1"
          "then we should stop tunting and choose a best alpha2 ever got";
  if (obj_val2 > obj_val1) {
    DEBUG() << "UNABLE TO TUNING A ALPHA2 BY HALFING";
    DEBUG() << "SELECT A BEST ALPHA2 EVER GOT";
    std::set<ParaPoint>::iterator it = tunedAlpha.begin();
    _alpha2 = it->first;
    _obj_val2 = it->second;

    _alpha3 = std::min(_alpha2 * 2, maxAlpha3);
    _obj_val3 = calobjval(grad, alpha3);


    DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
    DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;

    bool toParabolicFit = false;
    return toParabolicFit;
  }

  TRACE() << "now we can make sure that obj_val2 < obj_val1";

  const float alpha1 = 0;
  float linearFitAlph3 = (obj_val2 - obj_val1) / (alpha2 - alpha1) * (alpha3 - alpha1) + obj_val1;
  DEBUG() << __FUNCTION__ << format(" linear fit alpha3 = %e ") % linearFitAlph3;

  TRACE() << "keep the alpha we tuned";
  tunedAlpha.clear();
  tunedAlpha.insert(std::make_pair(alpha3, obj_val3));

  while (obj_val3 < linearFitAlph3 && obj_val3 < obj_val1 && alpha3 < maxAlpha3) {
    TRACE() << "if in this case, we should enlarge alpha3";
    alpha2 = alpha3;
    obj_val2 = obj_val3;

    alpha3 = std::min(alpha3 * 2, maxAlpha3);
    obj_val3 = calobjval(grad, alpha3);

    tunedAlpha.insert(std::make_pair(alpha3, obj_val3));

    DEBUG() << __FUNCTION__ << format(" tune alpha3, alpha2 = %e, obj_val2 = %e") % alpha2 % obj_val2;
    DEBUG() << __FUNCTION__ << format(" tune alpha3, alpha3 = %e, obj_val3 = %e") % alpha3 % obj_val3;
  }


  TRACE() << "If we couldnot tune a good alpha3";
  if (alpha3 > maxAlpha3 + 0.1) {
    DEBUG() << "UNABLE TO TUNING A ALPHA3 BY DOUBLING";
    DEBUG() << "SELECT A BEST ALPHA3 EVER GOT";
    std::set<ParaPoint>::iterator it = tunedAlpha.begin();
    _alpha3 = it->first;
    _obj_val3 = it->second;

    _alpha2 = _alpha3 / 2;
    _obj_val2 = calobjval(grad, alpha2);

    bool toParabolicFit = false;

    DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
    DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;
    return toParabolicFit;
  }

  /// return objval2 and objval3
  bool toParabolicFit = true;
  _alpha2 = alpha2;
  _alpha3 = alpha3;
  _obj_val2 = obj_val2;
  _obj_val3 = obj_val3;

  DEBUG() << __FUNCTION__ << format(" alpha2 = %e, obj_val2 = %e") % _alpha2 % _obj_val2;
  DEBUG() << __FUNCTION__ << format(" alpha3 = %e, obj_val3 = %e") % _alpha3 % _obj_val3;

  return toParabolicFit;
}

void UpdateSteplenOp::calsteplen(const std::vector<float>& grad,
    float obj_val1, int iter, float lambdaX, float lambdaZ, float &steplen, float &objval) {

    this->lambdaX = lambdaX;
    this->lambdaZ = lambdaZ;

  float dt = fmMethod.getdt();
  float dx = fmMethod.getdx();

  /// "calculate the initial value of alpha2 and alpha3";
  float max_alpha2, max_alpha3;

  calMaxAlpha2_3(fmMethod.getVelocity(), &grad[0], dt, dx, maxdv, max_alpha2, max_alpha3);
  DEBUG() << format("               max_alpha2 = %e,  max_alpha3: = %e") % max_alpha2 % max_alpha3;

  float alpha1 = 0, alpha2, alpha3;
  initAlpha23(max_alpha3, alpha2, alpha3);
  DEBUG() << format("after init alpha,  alpha2 = %e,      alpha3: = %e") % alpha2 % alpha3;

  float obj_val2, obj_val3;

  bool toParabolic = refineAlpha(grad, obj_val1, max_alpha3, alpha2, obj_val2, alpha3, obj_val3);

  float alpha4, obj_val4;
  if (toParabolic) {
    DEBUG() << "parabolic fit";

    TRACE() << "max_alpha3: " << max_alpha3;
    TRACE() << "alpha1: " << alpha1 << ", obj_val1: " << obj_val1 << ", alpha2: " << alpha2 << ", obj_val2: " << obj_val2
            << "alpha3: " << alpha3 << ", obj_val3: " << obj_val3;
    parabolaVertex(alpha1, obj_val1, alpha2, obj_val2, alpha3, obj_val3, max_alpha3, alpha4, obj_val4);
    TRACE() << "parabolaVertex done";

    if (alpha4 > max_alpha3) {
      DEBUG() << format("alpha4 = %e, max_alpha3 = %e") % alpha4 % max_alpha3;
      DEBUG() << format("alpha4 is greater than max_alpha3, set it to alpha3");
      alpha4 = max_alpha3;
    }
  } else {
    DEBUG() << "NO need to perform parabolic fit";
    alpha4 = alpha3;
    obj_val4 = obj_val3;
  }

  INFO() << format("iter %d  alpha  = %e total obj_val1 = %e") % iter % alpha1 % obj_val1;
  INFO() << format("iter %d  alpha2 = %e total obj_val2 = %e") % iter % alpha2 % obj_val2;
  INFO() << format("iter %d  alpha3 = %e total obj_val3 = %e") % iter % alpha3 % obj_val3;
  INFO() << format("iter %d  alpha4 = %e total obj_val4 = %e\n") % iter % alpha4 % obj_val4;

  preservedAlpha.alpha = alpha4;
  steplen = alpha4;
  objval = obj_val4;
}

void UpdateSteplenOp::bindEncSrcObs(const std::vector<float>& encsrc,
    const std::vector<float>& encobs) {
  this->encsrc = &encsrc;
  this->encobs = &encobs;
}

void UpdateSteplenOp::initAlpha23(float maxAlpha3, float &initAlpha2, float &initAlpha3) {
  const float minAlpha   = 1.0E-7;
  const float resetAlpha = 1.0E-4;

  if (!preservedAlpha.init) {
    preservedAlpha.init = true;
    preservedAlpha.alpha = maxAlpha3;
  }

  initAlpha3 = preservedAlpha.alpha;
  initAlpha3 = initAlpha3 < minAlpha ? resetAlpha : initAlpha3;
  initAlpha2 = initAlpha3 * 0.5;
}
