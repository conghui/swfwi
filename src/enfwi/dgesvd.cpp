/*
 * dgesvd.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

#include <algorithm>
#include "dgesvd.h"

extern "C" {
//#include <f2c.h>
//#include <clapack.h>
int dgesvd_(char *jobu, char *jobvt, int *m, int *n,
      double *a, int *lda, double *s, double *u, int *ldu, double *vt,
      int *ldvt, double *superb, int *lwork, int *info);
}

int LAPACKE_dgesvd_col_major(
    char jobu, char jobvt, int _m, int _n,
    double* a, int _lda, double* s, double* u, int _ldu, double* vt,
    int _ldvt, double* superb, int _lwork) {

  int m    = _m;
  int n    = _n;
  int lda  = _lda;
  int ldu  = _ldu;
  int ldvt = _ldvt;
  int lwork = _lwork;

  int info;
  dgesvd_(&jobu, &jobvt, &m, &n,
      a, &lda, s, u, &ldu, vt,
      &ldvt, superb, &lwork, &info);

  return info;
}
