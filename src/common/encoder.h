/*
 * encoder.h
 *
 *  Created on: Feb 27, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_ENCODER_H_
#define SRC_COMMON_ENCODER_H_

#include <vector>

class Encoder {
public:
  explicit Encoder(const std::vector<int> &code);
  std::vector<float> encodeSource(const std::vector<float> &wlt);
  std::vector<float> encodeObsData(const std::vector<float> &dobs, int nt, int ng);

private:
  const std::vector<int> &mCode;
};

#endif /* SRC_COMMON_ENCODER_H_ */
