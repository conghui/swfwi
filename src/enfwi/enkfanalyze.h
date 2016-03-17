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
#include "Matrix.h"

class EnkfAnalyze {
public:
  EnkfAnalyze(const Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, float sigmafactor);

  void analyze(std::vector<float *> &velSet) const;
  std::vector<float> createAMean(const std::vector<float *> &velSet) const;

protected:
  Matrix calGainMatrix(const std::vector<float *> &velSet) const;
  double initPerturbSigma(double maxHAP, float factor) const;
  void initGamma(const Matrix &perturbation, Matrix &gamma) const;
  void initPerturbation(Matrix &perturbation, double mean, double sigma) const;

protected:
  const Damp4t10d &fm;
  const std::vector<float> &wlt;
  const std::vector<float> &dobs;

  int modelSize;
  float sigmaFactor;
};

#endif /* SRC_ESS_FWI2D_ENKFANALYZE_H_ */
