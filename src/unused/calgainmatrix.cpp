/*
 * calgainmatrix.cpp
 *
 *  Created on: Mar 11, 2016
 *      Author: rice
 */

#include "calgainmatrix.h"

#include "logger.h"
#include <cassert>
#include <cmath>
#include <algorithm>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>
#include <boost/bind.hpp>
#include <string>
#include <cstdlib>
#include <vector>
#include <list>
#include <functional>
#include <cmath>
#include <numeric>
#include <cfloat>

#include <time.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "random-code.h"
#include "logger.h"
#include "Matrix.h"
#include "encoder.h"
#include "sfutil.h"
#include "common.h"
#include "dgesvd.h"

static const float enkf_sigma_factor = 0.2;

template <typename T>
float variance(const T *A_Begin, const T *A_End, const T *B_Begin) {
  float v = 0;

  for (; A_Begin != A_End; ++A_Begin, ++B_Begin) {
    v += (*A_Begin - *B_Begin) * (*A_Begin - *B_Begin);
  }

  return v;
}


void finalizeAMean(float *AMean) {
  free(AMean);
}

double initPerturbSigma(double maxHAP, float factor) {
  return maxHAP * factor;
}

void initPerturbation(Matrix &perturbation, double mean, double sigma) {
  //    int seed = time(NULL);
  int seed = 1;
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
  generator(boost::mt19937(seed),
            boost::normal_distribution<>(mean, sigma));

  std::generate(perturbation.getData(), perturbation.getData() + perturbation.size(), generator);
}



void initGamma(const Matrix &perturbation, Matrix &gamma) {
  assert(perturbation.isCompatible(gamma));

  const Matrix::value_type *p = perturbation.getData();
  for (int irow = 0; irow < gamma.getNumRow(); irow++) {
    Matrix::value_type sum = 0.0f;
    for (int icol = 0; icol < gamma.getNumCol(); icol++) {
      sum += p[icol * gamma.getNumRow() + irow];
    }

    Matrix::value_type avg = sum / gamma.getNumCol();

    for (int icol = 0; icol < gamma.getNumCol(); icol++) {
      int idx = icol * gamma.getNumRow() + irow;
      gamma.getData()[idx] = p[idx] - avg;
    }
  }
}

static void forwardModeling(const Damp4t10d &fmMethod,
    const std::vector<float> &encSrc,
    std::vector<float> &dobs /* output (fast: ng, slow: nt) */)
{
    int nx = fmMethod.getnx();
    int nz = fmMethod.getnz();
    int ns = fmMethod.getns();
    int ng = fmMethod.getng();
    int nt = fmMethod.getnt();

    std::vector<float> p0(nz * nx, 0);
    std::vector<float> p1(nz * nx, 0);

    for(int it=0; it<nt; it++) {

      fmMethod.addEncodedSource(&p1[0], &encSrc[it * ns]);

      fmMethod.stepForward(&p0[0], &p1[0]);

      std::swap(p1, p0);

      fmMethod.recordSeis(&dobs[it*ng], &p0[0]);


    }

}

Matrix calGainMatrix(const Damp4t10d &fm, const std::vector<float> &wlt, const std::vector<float> &dobs, std::vector<float *> &velSet,
                     std::vector<float> &resdSet, float *AMean, int modelSize, const int iter) {
  int N = velSet.size();
  int nt = fm.getnt();
  int ng = fm.getng();
//  int modelSize = fm.getnx() * fm.getnz();
  std::vector<std::vector<int> > codes(N);

  int numDataSamples = nt * ng;
  DEBUG() << "total samples in receivers: " << numDataSamples;
  Matrix HOnA(N, numDataSamples);
  Matrix HOnAMean(N, numDataSamples);
  Matrix D(N, numDataSamples);

  TRACE() << "use different codes for different column";
  int diffColDiffCodes = 0;
  int use_H_on_Amean = 0;
  INFO() << (diffColDiffCodes ? "different codes" : "same codes") << " for different velocity";
  INFO() << (use_H_on_Amean ? "" : "Not") << " apply forward modeling on AMean";

  for (int i = 0; i < N; i++) {
    int ns = fm.getns();
    DEBUG() << format("calculate HA on velocity %2d/%d") % (i + 1) % velSet.size();

    TRACE() << "making encoded shot, both sources and receivers";
    codes[i] = RandomCode::genPlus1Minus1(ns);
    std::copy(codes[i].begin(), codes[i].end(), std::ostream_iterator<int>(std::cout, ", ")); std::cout << "\n";

    Encoder encoder(codes[0]);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    {
      char buf[256];
      sprintf(buf, "encsrc%d.rsf", i);
//      sfFloatWrite1d(buf, &encsrc[0], encsrc.size());

      sprintf(buf, "encobs%d.rsf", i);
//      sfFloatWrite2d(buf, &encobs[0], ng, nt);
    }

    TRACE() << "save encoded data";
    std::vector<float> trans(encobs.size());
    matrix_transpose(&encobs[0], &trans[0], ng, nt);
    const float *pdata = &trans[0];
    std::copy(pdata, pdata + numDataSamples, D.getData() + i * numDataSamples);

    DEBUG() << format("sum D %.20f") % getSum(D);

    TRACE() << "H operate on A, and store data in HOnA";
    std::vector<float> dcal(encobs.size(), 0);

    ///// this decrease the performance ///
    Damp4t10d newfm = fm;
    Velocity curvel(std::vector<float>(velSet[i], velSet[i] + modelSize), fm.getnx(), fm.getnz());
    newfm.bindVelocity(curvel);
    forwardModeling(newfm, encsrc, dcal);
    matrix_transpose(&dcal[0], &trans[0], ng, nt);
    std::copy(pdata, pdata + numDataSamples, HOnA.getData() + i * numDataSamples);

    DEBUG() << format("sum HonA %.20f") % getSum(HOnA);

    TRACE() << "save the resd";
    const Matrix::value_type *obsDataBegin = D.getData() + i * numDataSamples;
    const Matrix::value_type *obsDataEnd   = obsDataBegin + numDataSamples;
    const Matrix::value_type *synDataBegin = HOnA.getData() + i * numDataSamples;
    resdSet[i] = variance(obsDataBegin, obsDataEnd, synDataBegin);

//    DEBUG() << format("resdSet[%d] %.20f") % i % resdSet[i];
  }

  DEBUG() << "sum of D: " << getSum(D);
  DEBUG() << "sum of HOnA: " << getSum(HOnA);
  DEBUG() << "sum of HOnAMean: " << getSum(HOnAMean);

  Matrix HA_Perturb(N, numDataSamples);
  initGamma(HOnA, HA_Perturb);
  DEBUG() << "sum of HA_Perturb: " << getSum(HA_Perturb);

  TRACE() << "initialize the perturbation";
  double mean = 0.0f;
  double maxHAP = std::abs(*std::max_element(HA_Perturb.getData(), HA_Perturb.getData() + HA_Perturb.size(), abs_less<float>));
  double sigma = initPerturbSigma(maxHAP, enkf_sigma_factor);
  DEBUG() << "sigma: " << sigma;
  Matrix perturbation(N, numDataSamples);
  initPerturbation(perturbation, mean, sigma);

  DEBUG() << "sum of perturbation: " << getSum(perturbation);

  TRACE() << "add perturbation to observed data";
  std::transform(D.getData(), D.getData() + D.size(), perturbation.getData(), D.getData(), std::plus<Matrix::value_type>());
  DEBUG() << "sum of D with perturbation added: " << getSum(D);

  TRACE() << "calculate the gamma";
  Matrix gamma(N, numDataSamples);
  initGamma(perturbation, gamma);
  DEBUG() << "sum of gamma: " << getSum(gamma);

  TRACE() << "calculate the band";
  Matrix band(N, numDataSamples);
  A_plus_B(HA_Perturb, gamma, band);
  DEBUG() << "sum of band: " << getSum(band);

  TRACE() << "svd of band";
  Matrix matU(N, numDataSamples);
  Matrix matS(1, N);
  Matrix matVt(N, N);
  Matrix superb(1, N - 1);

  int info = LAPACKE_dgesvd_col_major('S', 'S',
                                 band.getNumRow(), band.getNumCol(),
                                 band.getData(), band.getNumRow(),
                                 matS.getData(),
                                 matU.getData(), matU.getNumRow(),
                                 matVt.getData(), matVt.getNumRow(),
                                 superb.getData());

  /* Check for convergence */
  if ( info > 0 ) {
    ERROR() << "The algorithm computing SVD failed to converge.";
    exit( 1 );
  }

  DEBUG() << "sum of matU: " << getSum(matU);
  DEBUG() << "sum of matS: " << getSum(matS);

  TRACE() << "calculate matSSqInv";
  int clip = clipPosition(matS);
  DEBUG() << "clip: " << clip;
  DEBUG() << "matS: ";
  matS.print();

  Matrix matSSqInv(N, N);
  {
    Matrix::value_type *p = matSSqInv.getData();
    Matrix::value_type *s = matS.getData();
    const int nrow = matSSqInv.getNumRow();
    for (int i = 0; i < clip; i++) {
      p[i * nrow + i] = (1 / s[i]) * (1 / s[i]);
    }
  }
  DEBUG() << "sum of matSSqInv: " << getSum(matSSqInv);
  TRACE() << "matSSQInv: ";

  Matrix t0(N, numDataSamples); /// D - HA
  A_minus_B(D, HOnA, t0);
  DEBUG() << "sum of t0: " << getSum(t0);

  Matrix t1(N, N); /// U' * (D - HA)
  alpha_ATrans_B_plus_beta_C(1, matU, t0, 0, t1);
  DEBUG() << "sum of t1: " << getSum(t1);

  Matrix t2(N, N); /// HA' * U
  alpha_ATrans_B_plus_beta_C(1, HA_Perturb, matU, 0, t2);
  DEBUG() << "sum of t2: " << getSum(t2);

  Matrix t3(N, N); /// HA' * U * SSqInv
  alpha_A_B_plus_beta_C(1, t2, matSSqInv, 0, t3);
  DEBUG() << "sum of t3: " << getSum(t3);

  Matrix t4(N, N); /// HA' * U * SSqInv * U' * (D - HA)
  alpha_A_B_plus_beta_C(1, t3, t1, 0, t4);
  DEBUG() << "sum of t4: " << getSum(t4);
  DEBUG() << "print HA' * U * SSqInv * U' * (D - HA)";
  t4.print();

  finalizeAMean(AMean);
  return t4;

} /// end of function
