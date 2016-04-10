/*
 * dgesvd.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

#include <algorithm>
#include "dgesvd.h"

extern "C" {
#include <f2c.h>
#include <clapack.h>
}

int LAPACKE_dgesvd_col_major(
    char jobu, char jobvt, int _m, int _n,
    double* a, int _lda, double* s, double* u, int _ldu, double* vt,
    int _ldvt, double* superb, int _lwork) {

  integer m    = _m;
  integer n    = _n;
  integer lda  = _lda;
  integer ldu  = _ldu;
  integer ldvt = _ldvt;
  integer lwork = _lwork;

  integer info;
  dgesvd_(&jobu, &jobvt, &m, &n,
      a, &lda, s, u, &ldu, vt,
      &ldvt, superb, &lwork, &info);

  return info;
}

