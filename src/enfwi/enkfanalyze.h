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
#include "pMatrix.h"
#include <iostream>
#include "random-code.h"

class EnkfAnalyze {
public:
  EnkfAnalyze(const Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, float sigmafactor);

  void analyze(std::vector<float *> &totalVelSet, std::vector<float *> &velSet) const;
  void pAnalyze(std::vector<float *> &velSet, Matrix &lambdaSet, Matrix &ratioSet) const;
  std::vector<float> createAMean(const std::vector<float *> &velSet) const;
  std::vector<float> pCreateAMean(const std::vector<float *> &velSet, const int N) const;
	void check(std::vector<float> a, std::vector<float> b);
  void initLambdaSet(const std::vector<float*>& velSet, Matrix& lambdaSet, const Matrix& ratioSet) const;

protected:
  Matrix calGainMatrix(const std::vector<float *> &velSet, std::vector<int> code) const;
  Matrix pCalGainMatrix(const std::vector<float *> &velSet, std::vector<int> code, std::vector<float> &resdSet) const;
  double initPerturbSigma(double maxHAP, float factor) const;
  void initGamma(const Matrix &perturbation, Matrix &gamma) const;
  void pInitGamma(const Matrix &perturbation, Matrix &gamma, const int nSamples) const;
  void checkMatrix(const Matrix &a, const Matrix &b) const;
//  void initPerturbation(Matrix &perturbation, double mean, double sigma) const;
  void initPerturbation(Matrix& perturbation, const Matrix &HA_Perturb) const;
  void pInitPerturbation(Matrix& perturbation, const Matrix &HA_Perturb, const int rank, const int nSamples) const;
  void pInitPerturbation2(Matrix& perturbation, const Matrix &HA_Perturb, const int rank, const int nSamples) const;
  void pInitRatioPerturb(const Matrix &ratioSet, Matrix &ratioPerturb, int nsamples) const;


protected:
  static const int ENKF_SEED = 2;

protected:
  const Damp4t10d &fm;
  const std::vector<float> &wlt;
  const std::vector<float> &dobs;
  mutable RandomCodes enkfRandomCodes;

  int modelSize;
  float sigmaFactor;

  mutable boost::variate_generator<boost::mt19937, boost::normal_distribution<> > *generator;
  mutable float sigmaIter0;
  mutable bool initSigma;
};

#endif /* SRC_ESS_FWI2D_ENKFANALYZE_H_ */
