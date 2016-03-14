/*
 * essfwiframework.h
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_
#define SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_

#include "damp4t10d.h"
#include "updatevelop.h"

class EssFwiFramework {
public:
  EssFwiFramework(Damp4t10d &fmMethod, const UpdateVelOp &updateVelOp, const std::vector<float> &wlt, const std::vector<float> &dobs);

  void epoch(int iter);
  void writeVel(sf_file file) const;
  void bindEncSrcObs(const std::vector<float> &encsrc, const std::vector<float> &encobs);

private:
  float calobjval(const std::vector<float> &grad, float steplen) const;
  bool refineAlpha(const std::vector<float> &grad, float obj_val1, float maxAlpha3, float &_alpha2, float &_obj_val2, float &_alpha3, float &_obj_val3) const;
  float calsteplen(const std::vector<float> &grad, float obj_val1, int iter);
  void initAlpha23(float maxAlpha3, float &initAlpha2, float &initAlpha3);

private: /// from contructor, keep reference
  Damp4t10d &fmMethod;
  const UpdateVelOp &updateVelOp;
  const std::vector<float> &wlt;  /// wavelet
  const std::vector<float> &dobs; /// actual observed data (nt*ng*ns)


private: /// propagate from other construction
  int ns;
  int ng;
  int nt;
  int nx;
  int nz;
  float dx;
  float dt;

private:
  static std::vector<float> g0;               /// gradient in previous step
  static std::vector<float> updateDirection;


private:
  const std::vector<float> *encsrc;     /// encoded sources
  const std::vector<float> *encobs;     /// encoded observations

private:
  struct PreservedAlpha {
    float alpha;
    bool init;

    PreservedAlpha() : alpha(0), init(false) {}
  } preservedAlpha;
};

#endif /* SRC_ESS_FWI2D_ESSFWIFRAMEWORK_H_ */
