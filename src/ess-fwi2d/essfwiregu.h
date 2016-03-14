/*
 * essfwiregu.h
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ESSFWIREGU_H_
#define SRC_ESS_FWI2D_ESSFWIREGU_H_

#include "damp4t10d.h"

class EssFwiRegu {
public:
  EssFwiRegu(Damp4t10d &fmMethod, const std::vector<float> &wlt, const std::vector<float> &dobs);

  void epoch(int iter, int ivel, float lambdaX, float lambdaZ);
  void writeVel(sf_file file) const;

private: /// from contructor, keep reference
  Damp4t10d &fmMethod;
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
  static std::vector<float> g0;          /// gradient in previous step
  static std::vector<float> updateDirection;
};

#endif /* SRC_ESS_FWI2D_ESSFWIREGU_H_ */

