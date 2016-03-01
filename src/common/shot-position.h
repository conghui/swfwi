/*
 * shot-position.h
 *
 *  Created on: Mar 1, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_SHOT_POSITION_H_
#define SRC_COMMON_SHOT_POSITION_H_

#include <vector>

class ShotPosition {
public:
  ShotPosition(int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz);
  ShotPosition clipRange(int begin, int end) const;
  ShotPosition clip(int idx) const;
  int getx(int idx) const;
  int getz(int idx) const;

public:
  int ns;

private:
  std::vector<int> pos;
  int nz;
};

#endif /* SRC_COMMON_SHOT_POSITION_H_ */
