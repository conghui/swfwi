/*
 * calgainmatrix.h
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

#ifndef SRC_FWI_CALGAINMATRIX_H_
#define SRC_FWI_CALGAINMATRIX_H_

#include <vector>
#include "Matrix.h"
#include "damp4t10d.h"
Matrix calGainMatrix(const Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, std::vector<float *> &velSet, std::vector<float> &resdSet, float *AMean, int modelSize, const int iter);

void initGamma(const Matrix &perturbation, Matrix &gamma);
void finalizeAMean(float *AMean);

#endif /* SRC_FWI_CALGAINMATRIX_H_ */
