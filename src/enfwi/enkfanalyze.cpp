/*
 * enkfanalyze.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: rice
 */

#include <boost/bind.hpp>

#include "enkfanalyze.h"
#include "logger.h"
#include "sum.h"
#include "common.h"
#include "random-code.h"
#include "encoder.h"
#include "dgesvd.h"
#include "ReguFactor.h"
#include "aux.h"

namespace {
//std::vector<float> createAMean(const std::vector<float *> &velSet, int modelSize) {
//  std::vector<float> ret(modelSize);
//
//  for (int i = 0; i < modelSize; i++) {
//    float sum = 0.0;
//    for (size_t j = 0; j < velSet.size(); j++) {
//      sum += velSet[j][i];
//    }
//
//    ret[i] = sum / velSet.size();
//  }
//
//  return ret;
//}

void initAPerturb(Matrix &matAPerturb, const std::vector<float *> &velSet,
    const std::vector<float> &AMean, int modelSize) {
  assert((size_t)matAPerturb.getNumCol() == velSet.size()); // the matrix is row-majored
  assert(matAPerturb.getNumRow() == modelSize);

  for (int i = 0; i < matAPerturb.getNumCol(); i++) {
    Matrix::value_type *p = matAPerturb.getData() + (i * matAPerturb.getNumRow());
    std::transform(velSet[i], velSet[i] + modelSize, AMean.begin(), p, std::minus<Matrix::value_type>());
  }
}

} /// end of name space


EnkfAnalyze::EnkfAnalyze(const Damp4t10d &fm, const std::vector<float> &wlt,
    const std::vector<float> &dobs, float sigmafactor) :
  fm(fm), wlt(wlt), dobs(dobs), sigmaFactor(sigmafactor), sigmaIter0(0), initSigma(false)
{
  modelSize = fm.getnx() * fm.getnz();
}

void EnkfAnalyze::analyze(std::vector<float*>& velSet, Matrix &lambdaSet, Matrix &ratioSet) const {
  float dt = fm.getdt();
  float dx = fm.getdx();

  std::vector<float *> &A = velSet; /// velSet <==> matA
  int N = velSet.size();
  std::vector<float> AMean = createAMean(A);
  DEBUG() << "sum of AMean: " << sum(AMean);

  Matrix A_Perturb(N, modelSize);
  initAPerturb(A_Perturb, A, AMean, modelSize);
  DEBUG() << "sum of A_Perturb: " << getSum(A_Perturb);

  std::vector<float> resdSet(N);
  Matrix gainMatrix = calGainMatrix(velSet, resdSet);

  Matrix t5(N, modelSize);
  alpha_A_B_plus_beta_C(1, A_Perturb, gainMatrix, 0, t5);
  DEBUG() << "sum of t5: " << getSum(t5);

  TRACE() << "add the update back to velocity model";
  for (size_t i = 0; i < velSet.size(); i++) {
    float *vel = velSet[i];

    DEBUG() << format("before vel recovery, velset[%2d/%d], min: %f, max: %f") % i % velSet.size() %
            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

    TRACE() << "transform velocity to original";
    std::transform(vel, vel + modelSize, vel, boost::bind(velRecover<float>, _1, dx, dt));

    DEBUG() << format("vel recoverty,       velset[%2d/%d], min: %f, max: %f") % i % velSet.size() %
            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

    TRACE() << "add value calculated from ENKF to velocity";
    Matrix::value_type *pu = t5.getData() + i * t5.getNumRow();
    std::transform(vel, vel + modelSize, pu, vel, std::plus<Matrix::value_type>());

    TRACE() << "sum velset " << i << ": " << std::accumulate(vel, vel + modelSize, 0.0f);
    DEBUG() << format("after plus ENKF,     velset[%2d/%d], min: %f, max: %f\n") % i % velSet.size() %
            (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

    std::transform(vel, vel + modelSize, vel, boost::bind(velTrans<float>, _1, dx, dt));
  }

  TRACE() << "updating ratioset";
  Matrix ratio_Perturb(N, 2);
  initRatioPerturb(ratioSet, ratio_Perturb);

  Matrix t6(N, 2);
  alpha_A_B_plus_beta_C(1, ratio_Perturb, gainMatrix, 0, t6);

  A_plus_B(ratioSet, t6, ratioSet);

  int nx = fm.getnx();
  int nz = fm.getnz();
  float lambdaX = 0;
  float lambdaZ = 0;
  TRACE() << "updating lamdaset";
  for (int i = 0; i < N; i++) {
    ReguFactor fac(velSet[i], nx, nz, lambdaX, lambdaZ);
    Matrix::value_type *pLambda = lambdaSet.getData() + i * lambdaSet.getNumRow();
    Matrix::value_type *pRatio  = ratioSet.getData()  + i * ratioSet.getNumRow();
    pLambda[0] = pRatio[0] * resdSet[i] / fac.getWx2();
    pLambda[1] = pRatio[1] * resdSet[i] / fac.getWz2();
    DEBUG() << "PRATIO[0] " << pRatio[0];
    DEBUG() << "PRESD[i] " << resdSet[i];
    DEBUG() << "PLAMBDA[0] " << pLambda[0];
    DEBUG() << "";
  }
}


Matrix EnkfAnalyze::calGainMatrix(const std::vector<float*>& velSet, std::vector<float> &resdSet) const {
  int N = velSet.size();
  int nt = fm.getnt();
  int ng = fm.getng();
  int numDataSamples = nt * ng;

  std::vector<std::vector<int> > codes(N);
  Matrix HOnA(N, numDataSamples);
  Matrix D(N, numDataSamples);

  DEBUG() << "total samples in receivers: " << numDataSamples;
  TRACE() << "use different codes for different column";

  for (int i = 0; i < N; i++) {
    int ns = fm.getns();
    DEBUG() << format("calculate HA on velocity %2d/%d") % (i + 1) % velSet.size();

    /// "making encoded shot, both sources and receivers";
    codes[i] = RandomCode::genPlus1Minus1(ns);
    Encoder encoder(codes[0]);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    /// "save encoded data";
    std::vector<float> trans(encobs.size());
    matrix_transpose(&encobs[0], &trans[0], ng, nt);
    const float *pdata = &trans[0];
    std::copy(pdata, pdata + numDataSamples, D.getData() + i * numDataSamples);

    DEBUG() << format("sum D %.20f") % getSum(D);

    /// "H operate on A, and store data in HOnA";
    std::vector<float> dcal(encobs.size(), 0);

    ///// this decrease the performance ///
    Damp4t10d newfm = fm;
    Velocity curvel(std::vector<float>(velSet[i], velSet[i] + modelSize), fm.getnx(), fm.getnz());
    newfm.bindVelocity(curvel);

    DEBUG() << format("   curvel %.20f") % sum(curvel.dat);
    newfm.EssForwardModeling(encsrc, dcal);
    matrix_transpose(&dcal[0], &trans[0], ng, nt);
    std::copy(pdata, pdata + numDataSamples, HOnA.getData() + i * numDataSamples);

    DEBUG() << format("   sum HonA %.20f") % getSum(HOnA);

    TRACE() << "save the resd";
    const Matrix::value_type *obsDataBegin = D.getData() + i * numDataSamples;
    const Matrix::value_type *obsDataEnd   = obsDataBegin + numDataSamples;
    const Matrix::value_type *synDataBegin = HOnA.getData() + i * numDataSamples;
    resdSet[i] = variance(obsDataBegin, obsDataEnd, synDataBegin);
  }

  DEBUG() << "sum of D: " << getSum(D);
  DEBUG() << "sum of HOnA: " << getSum(HOnA);

  Matrix HA_Perturb(N, numDataSamples);
  initGamma(HOnA, HA_Perturb);
  DEBUG() << "sum of HA_Perturb: " << getSum(HA_Perturb);

  TRACE() << "initialize the perturbation";
  Matrix perturbation(N, numDataSamples);
  initPerturbation(perturbation, HA_Perturb);

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

  return t4;
}

double EnkfAnalyze::initPerturbSigma(double maxHAP, float factor) const {
  return maxHAP * factor;
}

void EnkfAnalyze::initGamma(const Matrix& perturbation, Matrix& gamma) const {
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

std::vector<float> EnkfAnalyze::createAMean(const std::vector<float*>& velSet) const {
  int modelSize = fm.getnx() * fm.getnz();
  std::vector<float> ret(modelSize);

  for (int i = 0; i < modelSize; i++) {
    float sum = 0.0;
    for (size_t j = 0; j < velSet.size(); j++) {
      sum += velSet[j][i];
    }

    ret[i] = sum / velSet.size();
  }

  return ret;
}

void EnkfAnalyze::initPerturbation(Matrix& perturbation, const Matrix &HA_Perturb) const {
  if (!initSigma) {
    initSigma = true;
    int seed = 1;
    double mean = 0;
    double maxHAP = std::abs(*std::max_element(HA_Perturb.getData(), HA_Perturb.getData() + HA_Perturb.size(), abs_less<float>));
    double sigma = initPerturbSigma(maxHAP, sigmaFactor);
    DEBUG() << "sigmaIter0: " << sigma;
    sigmaIter0 = sigma;
    generator = new boost::variate_generator<boost::mt19937, boost::normal_distribution<> >(boost::mt19937(seed), boost::normal_distribution<>(mean, sigma));
  }

  std::generate(perturbation.getData(), perturbation.getData() + perturbation.size(), *generator);

}

void EnkfAnalyze::initLambdaSet(const std::vector<float*>& velSet, Matrix& lambdaSet, const Matrix& ratioSet) const {
  DEBUG() << "initial lambda ratio is: " << ratioSet.getData()[0];

  int ng = fm.getng();
  int nt = fm.getnt();
  int ns = fm.getns();
  int nx = fm.getnx();
  int nz = fm.getnz();
  int modelSize = nx * nz;
  int diffColDiffCodes = 0;

  int N = velSet.size();
  std::vector<std::vector<int> > codes(N);

  int numDataSamples = ng * nt;
  std::vector<float> obsData(numDataSamples);
  std::vector<float> synData(numDataSamples);

  for (int i = 0; i < N; i++) {
    DEBUG() << format("init Lambda on velocity %2d/%d") % (i + 1) % velSet.size();

    ///  "making encoded shot, both sources and receivers";
    codes[i] = RandomCode::genPlus1Minus1(ns);
    Encoder encoder(codes[0]);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    TRACE() << "save encoded data";
    std::vector<float> trans(encobs.size());
    matrix_transpose(&encobs[0], &trans[0], ng, nt);
    const float *p= &trans[0];
    std::copy(p, p + numDataSamples, obsData.begin());

    std::vector<float> dcal(encobs.size(), 0);

    ///// this decrease the performance ///
    Damp4t10d newfm = fm;
    Velocity curvel(std::vector<float>(velSet[i], velSet[i] + modelSize), fm.getnx(), fm.getnz());
    newfm.bindVelocity(curvel);
    newfm.EssForwardModeling(encsrc, dcal);
    matrix_transpose(&dcal[0], &trans[0], ng, nt);
    std::copy(p, p + numDataSamples, synData.begin());

    TRACE() << "calculate the data residule";
    float resd = variance(obsData, synData);

    TRACE() << "create regularization factor";
    ReguFactor reguFac(velSet[i], nx, nz);
    float Wx2 = reguFac.getWx2();
    float Wz2 = reguFac.getWz2();

    TRACE() << "update LambdaSet";
    assert(lambdaSet.getNumCol() == N);
    Matrix::value_type *plamda = lambdaSet.getData() + i * lambdaSet.getNumRow();
    const Matrix::value_type *pratio = ratioSet.getData() + i * ratioSet.getNumRow();
    plamda[0] = pratio[0] * resd / Wx2; /// lambdaX
    plamda[1] = pratio[1] * resd / Wz2; /// lamdaZ
    DEBUG() << "pratio[0] " << pratio[0];
    DEBUG() << "plambda[0] " << plamda[0];

  }
}

void EnkfAnalyze::initRatioPerturb(const Matrix& ratioSet,
    Matrix& ratioPerturb) const {
  initGamma(ratioSet, ratioPerturb);
}
