/*
 * sfutil.h
 *
 *  Created on: Mar 3, 2016
 *      Author: rice
 */



#ifndef SRC_UTIL_SFUTIL_H_
#define SRC_UTIL_SFUTIL_H_

void sfFloatWrite1d(const char *fn, const float *dat, int n1, float d1 = 1, float o1 = 0);
void sfFloatWrite2d(const char *fn, const float *dat, int n1, int n2,
    float d1 = 1, float d2 = 1, float o1 = 0, float o2 = 0);

void sfDoubleWrite2d(const char *fn, const double *dat, int n1, int n2,
    float d1 = 1, float d2 = 1, float o1 = 0, float o2 = 0);
#endif /* SRC_UTIL_SFUTIL_H_ */
