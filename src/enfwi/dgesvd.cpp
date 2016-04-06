/*
 * dgesvd.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

//#include <mkl.h>
#include <algorithm>
#include "dgesvd.h"
extern "C" {
#include <rsf.h>
}

/**
 *     int info = LAPACKE_dgesvd_col_major('S', 'S',
        band.getNumRow(), band.getNumCol(),
        band.getData(), band.getNumRow(),
        matS.getData(),
        matU.getData(), matU.getNumRow(),
        matVt.getData(), matVt.getNumRow(),
        superb.getData())
 */
int LAPACKE_dgesvd_col_major(char jobu, char jobvt, int m,
    int n, double* a, int lda, double* s, double* u, int ldu, double* vt,
    int ldvt, double* superb, int lwork) {

//  return LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, m,
//    n, a, lda, s, u, ldu, vt,
//    ldvt, superb);

//  int lwork = std::max(1, std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n)));
//  int lwork =

  int info;
  dgesvd_(&jobu, &jobvt, &m,
      &n, a, &lda, s, u, &ldu, vt,
      &ldvt, superb, &lwork, &info);

  return info;

//           dgesvd_(&jobu, &jobvt,
//                                         &numDataSamples, &N,
//                                         band, &numDataSamples,
//                                         matS,
//                                         matU, &numDataSamples,
//                                         matVt, &N,
//                                         work, &lwork, &info);
}

