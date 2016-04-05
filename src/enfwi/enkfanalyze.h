/*
 * enkfanalyze.h
 *
 *  Created on: Mar 15, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_ENKFANALYZE_H_
#define SRC_ESS_FWI2D_ENKFANALYZE_H_

#include <vector>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include "damp4t10d.h"
#include "Matrix.h"

class EnkfAnalyze {
public:
  EnkfAnalyze(const Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, float sigmafactor);

  void analyze(std::vector<float *> &totalVelSet, std::vector<float *> &velSet) const;
  std::vector<float> createAMean(const std::vector<float *> &velSet) const;

protected:
  Matrix calGainMatrix(const std::vector<float *> &velSet) const;
  double initPerturbSigma(double maxHAP, float factor) const;
  void initGamma(const Matrix &perturbation, Matrix &gamma) const;
//  void initPerturbation(Matrix &perturbation, double mean, double sigma) const;
  void initPerturbation(Matrix& perturbation, const Matrix &HA_Perturb) const;

protected:
  const Damp4t10d &fm;
  const std::vector<float> &wlt;
  const std::vector<float> &dobs;

  int modelSize;
  float sigmaFactor;

  mutable boost::variate_generator<boost::mt19937, boost::normal_distribution<> > *generator;
  mutable float sigmaIter0;
  mutable bool initSigma;
};

#endif /* SRC_ESS_FWI2D_ENKFANALYZE_H_ */
