/*
 * dgesvd.h
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

#ifndef SRC_FWI_DGESVD_H_
#define SRC_FWI_DGESVD_H_


int LAPACKE_dgesvd_col_major(char jobu, char jobvt,
   int m, int n, double* a,
   int lda, double* s, double* u, int ldu,
   double* vt, int ldvt, double* superb, int lwork);

#endif /* SRC_FWI_DGESVD_H_ */
