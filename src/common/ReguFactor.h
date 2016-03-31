/*
 * ReguFactor.h
 *
 *  Created on: Mar 4, 2015
 *      Author: Conghui He (heconghui@gmail.com)
 */

#ifndef REGUFACTOR_H_
#define REGUFACTOR_H_

#include <vector>

class ReguFactor {
 public:
  ReguFactor(const float *vel, int nx, int nz, float lambdaX = 0, float lambdaZ = 0);

  float getWx2() const;
  float getWz2() const;

  float getReguTerm();
  const float *getReguGradient();

 private:
  const float *mModel;
  int mNx;
  int mNz;

  float mLambdaX;
  float mLambdaZ;

  float mReguTermValue;
  std::vector<float> mReguGradVector;
};

#endif /* REGUFACTOR_H_ */
