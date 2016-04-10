#include <cmath>
#include <vector>
#include <iostream>

#include <cmath>
#include "logger.h"

static void calParabolaVertex(float x1, float y1, float x2, float y2, float x3, float y3, float &xv, float &yv) {
  DEBUG() << "inside function: " << __PRETTY_FUNCTION__;
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

  DEBUG() << "after calParabolaVertex";

  //if (std::abs(k2 - k1) < 0.001 * (std::max(std::abs(k2), std::abs(k1))) ||
      //(xv == -std::numeric_limits<double>::quiet_NaN())) {
    //WARNING() << "THE SET OF POINTS DON'T FIT PARABOLIC WELL, SET y TO y3 ON PURPOSE JUST FOR INDICATION";
////    xv = std::min(2 * x3, max_alpha3);
////    yv = -std::numeric_limits<double>::quiet_NaN(); /// indicating what's happening
    //xv = x3;
    //yv = y3;
  //}
}

int main(int argc, char* argv[]) {

  std::cout << "hello, how are you?" << std::endl;
  for (int i = 0; i < argc; i++){
    printf("%d %s\n", i, argv[i]);
  }

  float max_alpha3 = 0.00012343;
  float alpha1 = 0;
  float alpha2 = 6.17151e-05;
  float alpha3 = 0.00012343;
  float alpha4;
  float obj_val1 = 6.84734e+06;
  float obj_val2 = 5.49986e+06;
  float obj_val3 = 5.26678e+06;
  float obj_val4;

  TRACE() << "alpha1: " << alpha1 << ", obj_val1: " << obj_val1 << ", alpha2: " << alpha2 << ", obj_val2: " << obj_val2
          << ", alpha3: " << alpha3 << ", obj_val3: " << obj_val3;
  parabolaVertex(alpha1, obj_val1, alpha2, obj_val2, alpha3, obj_val3, max_alpha3, alpha4, obj_val4);

  INFO() << "program exit normally";
  return 0;
}

