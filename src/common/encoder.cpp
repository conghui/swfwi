/*
 * encoder.cpp
 *
 *  Created on: Feb 27, 2016
 *      Author: rice
 */

#include "encoder.h"


Encoder::Encoder(const std::vector<int>& code) : mCode(code) {
}

std::vector<float> Encoder::encodeSource(const std::vector<float>& wlt) {
  int ns = mCode.size();
  int nt = wlt.size();
  std::vector<float> encSrc(ns * nt);

  for (int it = 0; it < nt; it++) {
    for (int is = 0; is < ns; is++) {
      encSrc[it * ns + is] = wlt[it] * mCode[is];
    }
  }

  return encSrc;
}

std::vector<float> Encoder::encodeObsData(const std::vector<float>& dobs,
    int nt, int ng) {
  std::vector<float> encObs(nt * ng, 0);
  int ns = mCode.size();

  for (int is = 0; is < ns; is++) {
    int oneShotDataSize = nt * ng;

    for (int j = 0; j < oneShotDataSize; j++) {
      encObs[j] += dobs[is * oneShotDataSize + j] * mCode[is];
    }
  }

  return encObs;
}
