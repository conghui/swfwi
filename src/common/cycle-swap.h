/*
 * cycle-swap.h
 *
 *  Created on: Feb 24, 2016
 *      Author: Conghui He
 */

#ifndef SRC_UTIL_CYCLE_SWAP_H_
#define SRC_UTIL_CYCLE_SWAP_H_

#include <algorithm>

template <typename T>
void cycleSwap(T &a, T &b, T &c) {
  std::swap(a, b);
  std::swap(b, c);
}



#endif /* SRC_UTIL_CYCLE_SWAP_H_ */
