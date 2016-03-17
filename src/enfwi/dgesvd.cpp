/*
 * dgesvd.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

#include <mkl.h>
#include "dgesvd.h"

int LAPACKE_dgesvd_col_major(char jobu, char jobvt, int m,
    int n, double* a, int lda, double* s, double* u, int ldu, double* vt,
    int ldvt, double* superb) {

  return LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m,
    n, a, lda, s, u, ldu, vt,
    ldvt, superb);
}
