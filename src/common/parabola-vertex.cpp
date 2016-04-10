/*
 * parabola-vertex.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#include <cmath>
#include "parabola-vertex.h"
#include "logger.h"

static void calParabolaVertex(float x1, float y1, float x2, float y2, float x3, float y3, float &xv, float &yv) {
  double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
  double A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
  double B     = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
  double C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

  xv = -B / (2 * A);
  yv = C - B * B / (4 * A);
}

void parabolaVertex(float x1, float y1, float x2, float y2, float x3, float y3, float max_alpha3, float &xv, float &yv) {
  double k2 = (y3 - y2) / (x3 - x2);
  double k1 = (y2 - y1) / (x2 - x1);

  calParabolaVertex(x1, y1, x2, y2, x3, y3, xv, yv);

  if(std::abs(k2-k1)<0.001*(std::max(std::abs(k2),std::abs(k1)))|| std::isnan(-xv)) {
    WARNING() << "THE SET OF POINTS DON'T FIT PARABOLIC WELL, SET y TO y3 ON PURPOSE JUST FOR INDICATION";
    xv = x3;
    yv = y3;
  }
}
