/*
 * random-code.cpp
 *
 *  Created on: Feb 27, 2016
 *      Author: rice
 */

#include <cstdlib>
#include "random-code.h"

std::vector<int> RandomCode::genPlus1Minus1(int nshots) {
  std::vector<int> codes(nshots);
  for (int i = 0; i < nshots; i++) {
    int v = (static_cast<float>(rand()) / RAND_MAX) > 0.5 ? 1 : -1;
    codes[i] = v;
  }
  return codes;
}
