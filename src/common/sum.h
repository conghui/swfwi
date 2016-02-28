/*
 * sum.h
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_SUM_H_
#define SRC_COMMON_SUM_H_

#include <vector>
#include <numeric>

template <typename T>
T sum(const std::vector<T> &v) {
  return std::accumulate(v.begin(), v.end(), static_cast<T>(0));
}

template <typename T>
T sum(const T *v, int size) {
  return std::accumulate(v, v + size, static_cast<T>(0));
}

#endif /* SRC_COMMON_SUM_H_ */
