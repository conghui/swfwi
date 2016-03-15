/*
 * enkfanalyze.h
 *
 *  Created on: Mar 15, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ENKFANALYZE_H_
#define SRC_ESS_FWI2D_ENKFANALYZE_H_

#include <vector>

#include "damp4t10d.h"

class EnkfAnalyze {
public:
  EnkfAnalyze();

  void analyze(std::vector<float *> &velSet);

private:
  const Damp4t10d &fm;

  int modelSize;
};

#endif /* SRC_ESS_FWI2D_ENKFANALYZE_H_ */
