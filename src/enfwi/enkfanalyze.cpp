/*
 * enkfanalyze.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: rice
 */

#include <boost/bind.hpp>
#include <mpi.h>

#include "enkfanalyze.h"
#include "logger.h"
#include "sum.h"
#include "common.h"
#include "random-code.h"
#include "encoder.h"
#include "dgesvd.h"

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

void EnkfAnalyze::analyze(std::vector<float*>& totalVelSet, std::vector<float *> &velSet) const {
  Matrix gainMatrix = calGainMatrix(velSet);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    float dt = fm.getdt();
    float dx = fm.getdx();

    std::vector<float *> &A = totalVelSet; /// velSet <==> matA
    int N = totalVelSet.size();
    std::vector<float> AMean = createAMean(A);
    DEBUG() << "sum of AMean: " << sum(AMean);

    Matrix A_Perturb(N, modelSize);
    initAPerturb(A_Perturb, A, AMean, modelSize);
    DEBUG() << "sum of A_Perturb: " << getSum(A_Perturb);

    Matrix t5(N, modelSize);
    alpha_A_B_plus_beta_C(1, A_Perturb, gainMatrix, 0, t5);
    DEBUG() << "sum of t5: " << getSum(t5);

    TRACE() << "add the update back to velocity model";
    for (size_t i = 0; i < totalVelSet.size(); i++) {
      float *vel = totalVelSet[i];

      DEBUG() << format("before vel recovery, velset[%2d/%d], min: %f, max: %f") % i % totalVelSet.size() %
          (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

      TRACE() << "transform velocity to original";
      std::transform(vel, vel + modelSize, vel, boost::bind(velRecover<float>, _1, dx, dt));

      DEBUG() << format("vel recoverty,       velset[%2d/%d], min: %f, max: %f") % i % totalVelSet.size() %
          (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

      TRACE() << "add value calculated from ENKF to velocity";
      Matrix::value_type *pu = t5.getData() + i * t5.getNumRow();
      std::transform(vel, vel + modelSize, pu, vel, std::plus<Matrix::value_type>());

      TRACE() << "sum velset " << i << ": " << std::accumulate(vel, vel + modelSize, 0.0f);
      DEBUG() << format("after plus ENKF,     velset[%2d/%d], min: %f, max: %f\n") % i % totalVelSet.size() %
          (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));

      std::transform(vel, vel + modelSize, vel, boost::bind(velTrans<float>, _1, dx, dt));
    }

  }
}

Matrix EnkfAnalyze::calGainMatrix(const std::vector<float*>& velSet) const {
  int local_n = velSet.size();
  int nt = fm.getnt();
  int ng = fm.getng();
  int numDataSamples = nt * ng;

  int N = 0;
  MPI_Reduce(&local_n, &N, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); /// reduce # of total samples to rank 0


  std::vector<int> code = RandomCode::genPlus1Minus1(fm.getns());
  MPI_Bcast(&code[0], code.size(), MPI_INT, 0, MPI_COMM_WORLD);   /// broadcast the code to all other processes

  Matrix local_HOnA(local_n, numDataSamples);
  Matrix local_D(local_n, numDataSamples);
  Matrix HOnA(N, numDataSamples);
  Matrix D(N, numDataSamples);

  for (int i = 0; i < local_n; i++) {
    DEBUG() << format("calculate HA on velocity %2d/%d") % (i + 1) % velSet.size();

    /// "making encoded shot, both sources and receivers";
    Encoder encoder(code);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    /// "save encoded data";
    std::vector<float> trans(encobs.size());
    matrix_transpose(&encobs[0], &trans[0], ng, nt);
    const float *pdata = &trans[0];
    std::copy(pdata, pdata + numDataSamples, local_D.getData() + i * numDataSamples);

    DEBUG() << format("sum D %.20f") % getSum(local_D);

    /// "H operate on A, and store data in HOnA";
    std::vector<float> dcal(encobs.size(), 0);

    ///// this decrease the performance ///
    Damp4t10d newfm = fm;
    Velocity curvel(std::vector<float>(velSet[i], velSet[i] + modelSize), fm.getnx(), fm.getnz());
    newfm.bindVelocity(curvel);

    DEBUG() << format("   curvel %.20f") % sum(curvel.dat);
    newfm.EssForwardModeling(encsrc, dcal);
    matrix_transpose(&dcal[0], &trans[0], ng, nt);
    std::copy(pdata, pdata + numDataSamples, local_HOnA.getData() + i * numDataSamples);

    DEBUG() << format("   sum HonA %.20f") % getSum(local_HOnA);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int count = local_D.getNumCol() * local_D.getNumRow();
  MPI_Gather(local_D.getData(), count, MPI_DOUBLE, D.getData(), count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(local_HOnA.getData(), count, MPI_DOUBLE, HOnA.getData(), count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  Matrix t4(N, N); /// HA' * U * SSqInv * U' * (D - HA)
  if (rank == 0) {
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

    int lwork = std::max(1, std::max(3 * std::min(numDataSamples, N) + std::max(numDataSamples, N), 5 * std::min(numDataSamples, N)));
    Matrix superb(1, lwork);

    int info = LAPACKE_dgesvd_col_major('S', 'S',
        band.getNumRow(), band.getNumCol(),
        band.getData(), band.getNumRow(),
        matS.getData(),
        matU.getData(), matU.getNumRow(),
        matVt.getData(), matVt.getNumRow(),
        superb.getData(), lwork);

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

    alpha_A_B_plus_beta_C(1, t3, t1, 0, t4);
    DEBUG() << "sum of t4: " << getSum(t4);
    DEBUG() << "print HA' * U * SSqInv * U' * (D - HA)";
    t4.print();
  }

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
