/*
 * preserved-alpha.h
 *
 *  Created on: Mar 7, 2016
 *      Author: rice
 */

#ifndef SRC_ESS_FWI2D_PRESERVED_ALPHA_H_
#define SRC_ESS_FWI2D_PRESERVED_ALPHA_H_

#include <vector>
#include <algorithm>

class PreservedAlpha {
 public:
  const static int MAX = 1024;

 public:
  static PreservedAlpha &instance() {
    static PreservedAlpha ins;
    return ins;
  }

  std::vector<float> &getAlpha() {
    return mAlpha;
  }

  std::vector<bool> &getIsInit() {
    return mIsInit;
  }

 private:
  PreservedAlpha() : mObjval(MAX, 1E11), mAlpha(MAX, 0), mIsInit(MAX, false) {
    std::fill(mObjval.begin(), mObjval.end(), 1E11);
  }

  ~PreservedAlpha() {}
  PreservedAlpha(const PreservedAlpha &);
  void operator=(const PreservedAlpha &);

 private:
  std::vector<float> mObjval;
  std::vector<float> mAlpha;
  std::vector<bool>  mIsInit;
};


#endif /* SRC_ESS_FWI2D_PRESERVED_ALPHA_H_ */
