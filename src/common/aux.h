/*
 *
 *  Created on: Nov 13, 2014
 *      Author: heconghui@gmail.com
 */

#ifndef STRUTIL_H_
#define STRUTIL_H_

#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cassert>
#include <iomanip>

template <typename T>
std::string to_str(T val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}

template <typename T>
std::string to_str(T v, int ndigit) {
  std::stringstream ss;

  std::fixed(ss);
  ss << std::setprecision(ndigit) << v;

  return ss.str();
}

template <typename T>
bool is_close(T a, T b, double epsilon = 1e-5) {
  return std::fabs(a - b) < epsilon;
}

template <typename T>
float variance(const std::vector<T> &A, const std::vector<T> &B) {
  float v = 0;

  assert(A.size() == B.size());
  for (size_t i = 0; i < A.size(); i++) {
    v += (A[i] - B[i]) * (A[i] - B[i]);
  }

  return v;
}

template <typename T>
void writeBin(const std::string &fn, T *data, int bytes) {
  std::ofstream ofs(fn.c_str());
  assert(ofs.good());
  ofs.write(reinterpret_cast<char *>(data), bytes);
  ofs.close();
}

template <typename T>
void readBin(const std::string &fn, T *data, int bytes) {
  std::ifstream ifs(fn.c_str());
  assert(ifs.good());
  ifs.read(reinterpret_cast<char *>(data), bytes);
  ifs.close();
}

template <typename T>
T addSquare(T x, T y) {
  return x + y * y;
}

//template <typename T>
//bool abs_less(T a, T b) {
//  return std::abs(a) < std::abs(b);
//}

#endif /* STRUTIL_H_ */
