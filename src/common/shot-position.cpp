/*
 * shot-position.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: rice
 */

#include "shot-position.h"
#include <algorithm>

ShotPosition::ShotPosition(int szbeg, int sxbeg, int jsz, int jsx, int _ns, int _nz) :
    ns(_ns), pos(_ns), nz(_nz)

{
  for (int is = 0; is < ns; is++) {
    int sz = szbeg + is * jsz;
    int sx = sxbeg + is * jsx;
    pos[is] = sz + nz * sx;
  }
}

int ShotPosition::getx(int idx) const {
  return pos[idx] / nz;
}

ShotPosition ShotPosition::clipRange(int begin, int end) const {
  ShotPosition ret = *this;
  ret.ns = end - begin + 1;
  ret.pos.resize(ret.ns);
  std::copy(&pos[begin], &pos[end + 1], &ret.pos[0]);

  return ret;
}

int ShotPosition::getz(int idx) const {
  return pos[idx] % nz;
}

ShotPosition ShotPosition::clip(int idx) const {
  return clipRange(idx, idx);
}
