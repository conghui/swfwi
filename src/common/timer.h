/*
 * timer.h
 *
 *  Created on: Mar 30, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_TIMER_H_
#define SRC_COMMON_TIMER_H_

#include <iostream>
#include <ctime>

class Timer
{
public:
  Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

  double elapsed() {
    clock_gettime(CLOCK_REALTIME, &end_);
    return end_.tv_sec - beg_.tv_sec +
        (end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
  }

  void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

private:
  timespec beg_, end_;
};



#endif /* SRC_COMMON_TIMER_H_ */
