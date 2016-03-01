/*
 * shot-position.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: rice
 */

#include "shot-position.h"


ShotPosition::ShotPosition(int szbeg, int sxbeg, int jsz, int jsx, int _ns, int _nz) :
    pos(_ns), ns(_ns), nz(_nz)

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

int ShotPosition::getz(int idx) const {
  return pos[idx] % nz;
}
