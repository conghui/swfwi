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
#include "aux.h"
#include "ReguFactor.h"

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

	//printf("initAPerturb, col = %d, row = %d\n", matAPerturb.getNumCol(), matAPerturb.getNumRow()); 
  for (int i = 0; i < matAPerturb.getNumCol(); i++) {
    Matrix::value_type *p = matAPerturb.getData() + (i * matAPerturb.getNumRow());
    std::transform(velSet[i], velSet[i] + modelSize, AMean.begin(), p, std::minus<Matrix::value_type>());
  }
}

} /// end of name space


EnkfAnalyze::EnkfAnalyze(const Damp4t10d &fm, const std::vector<float> &wlt,
    const std::vector<float> &dobs, float sigmafactor) :
  fm(fm), wlt(wlt), dobs(dobs), enkfRandomCodes(ENKF_SEED), sigmaFactor(sigmafactor), sigmaIter0(0), initSigma(false)
{
  modelSize = fm.getnx() * fm.getnz();
}

void EnkfAnalyze::analyze(std::vector<float*>& totalVelSet, std::vector<float *> &velSet) const {
  std::vector<int> code = enkfRandomCodes.genPlus1Minus1(fm.getns());
  Matrix gainMatrix = calGainMatrix(velSet, code);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0)
	{
		DEBUG() << "sum of gainMatrix: " << getSum(gainMatrix);
	}


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

void EnkfAnalyze::pAnalyze(std::vector<float *> &velSet, Matrix &lambdaSet, Matrix &ratioSet) const {
  std::vector<int> code = enkfRandomCodes.genPlus1Minus1(fm.getns());

  int local_n = velSet.size();
  std::vector<float> resdSet(local_n);
  Matrix pGainMatrix = pCalGainMatrix(velSet, code, resdSet);


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*
	gainMatrix.print("gainMatrix");
	char fname[20];
	sprintf(fname, "pGainMatrix%d", rank);
	pGainMatrix.print(fname);
   */


  int N = 0;
  MPI_Reduce(&local_n, &N, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); /// reduce # of total samples to rank 0
  int nSamples = N;
  MPI_Bcast(&nSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);

  Matrix::value_type sum_pGainMatrix = pGetSum(pGainMatrix, nSamples);
  if(rank == 0)
  {
    //DEBUG() << "sum of gainMatrix: " << getSum(gainMatrix);
    DEBUG() << "sum of pGainMatrix: " << sum_pGainMatrix;
  }


  float dt = fm.getdt();
  float dx = fm.getdx();
  std::vector<float *> &local_A = velSet; /// velSet <==> matA
  std::vector<float> pAMean = pCreateAMean(local_A, nSamples);
  Matrix local_A_Perturb(local_n, modelSize);
  initAPerturb(local_A_Perturb, local_A, pAMean, modelSize);
  /*
	char filename[20];
	sprintf(filename, "HA_Perturb%d.txt", rank);
	local_A_Perturb.print(filename);
   */
  Matrix::value_type sum_A_Perturb = pGetSum(local_A_Perturb, nSamples);
  Matrix local_t5(local_n, modelSize);
  pAlpha_A_B_plus_beta_C(1.0, local_A_Perturb, 1, pGainMatrix, 1, 0.0, local_t5, 1, nSamples);
  Matrix::value_type sum_local_t5 = pGetSum(local_t5, nSamples);

  if(rank == 0)
  {
    DEBUG() << "sum of pAMean: " << sum(pAMean);
    DEBUG() << "sum of local_A_Perturb: " << sum_A_Perturb;
    DEBUG() << "sum of pt5: " << sum_local_t5;
    TRACE() << "add the update back to velocity model";
  }
  for (size_t i = 0; i < velSet.size(); i++) {
    float *vel = velSet[i];
    int absvel = rank * local_n + i + 1;
    DEBUG() << format("before vel recovery, velset[%2d/%d], min: %f, max: %f") % absvel % nSamples %
        (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));
    TRACE() << "transform velocity to original";
    std::transform(vel, vel + modelSize, vel, boost::bind(velRecover<float>, _1, dx, dt));
    DEBUG() << format("vel recoverty,       velset[%2d/%d], min: %f, max: %f") % absvel % nSamples %
        (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));
    TRACE() << "add value calculated from ENKF to velocity";
    Matrix::value_type *pu = local_t5.getData() + i * local_t5.getNumRow();
    std::transform(vel, vel + modelSize, pu, vel, std::plus<Matrix::value_type>());
    TRACE() << "sum velset " << i << ": " << std::accumulate(vel, vel + modelSize, 0.0f);
    DEBUG() << format("after plus ENKF,     velset[%2d/%d], min: %f, max: %f\n") % absvel % nSamples %
        (*std::min_element(vel, vel + modelSize)) % (*std::max_element(vel, vel + modelSize));
    std::transform(vel, vel + modelSize, vel, boost::bind(velTrans<float>, _1, dx, dt));
  }

  TRACE() << "updating ratioset";
  Matrix ratio_Perturb(local_n, 2);
  pInitRatioPerturb(ratioSet, ratio_Perturb, nSamples);

  Matrix local_t6(local_n, 2);
  pAlpha_A_B_plus_beta_C(1, ratio_Perturb, 1, pGainMatrix, 1, 0, local_t6, 1, nSamples);

  A_plus_B(ratioSet, local_t6, ratioSet);

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

Matrix EnkfAnalyze::calGainMatrix(const std::vector<float*>& velSet, std::vector<int> code) const {
  int local_n = velSet.size();
  int nt = fm.getnt();
  int ng = fm.getng();
  int numDataSamples = nt * ng;

  int N = 0;
  MPI_Reduce(&local_n, &N, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); /// reduce # of total samples to rank 0

	if(code.size() == 0)
	{
		code = enkfRandomCodes.genPlus1Minus1(fm.getns());
	}
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

		//band.print("band.txt");

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

Matrix EnkfAnalyze::pCalGainMatrix(const std::vector<float*>& velSet, std::vector<int> code, std::vector<float> &resdSet) const {
  int local_n = velSet.size();
  int nt = fm.getnt();
  int ng = fm.getng();
  int numDataSamples = nt * ng;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int N = 0;
  MPI_Reduce(&local_n, &N, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); /// reduce # of total samples to rank 0

	int nSamples = N;
	MPI_Bcast(&nSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(code.size() == 0)
	{
		std::vector<int> code = enkfRandomCodes.genPlus1Minus1(fm.getns());
	}
  MPI_Bcast(&code[0], code.size(), MPI_INT, 0, MPI_COMM_WORLD);   /// broadcast the code to all other processes

  Matrix local_HOnA(local_n, numDataSamples);
  Matrix local_D(local_n, numDataSamples);

  for (int i = 0; i < local_n; i++) {
		int absvel = rank * local_n + i + 1;
    DEBUG() << format("calculate HA on velocity %2d/%d") % absvel % nSamples;

    /// "making encoded shot, both sources and receivers";
    Encoder encoder(code);
    std::vector<float> encobs = encoder.encodeObsData(dobs, nt, ng);
    std::vector<float> encsrc  = encoder.encodeSource(wlt);

    /// "save encoded data";
    std::vector<float> trans(encobs.size());
    matrix_transpose(&encobs[0], &trans[0], ng, nt);
    const float *pdata = &trans[0];
    std::copy(pdata, pdata + numDataSamples, local_D.getData() + i * numDataSamples);

    DEBUG() << format("parallel: sum D %.20f") % getSum(local_D);

    /// "H operate on A, and store data in HOnA";
    std::vector<float> dcal(encobs.size(), 0);

    ///// this decrease the performance ///
    Damp4t10d newfm = fm;
    Velocity curvel(std::vector<float>(velSet[i], velSet[i] + modelSize), fm.getnx(), fm.getnz());
    newfm.bindVelocity(curvel);

    DEBUG() << format("parallel: curvel %.20f") % sum(curvel.dat);
    newfm.EssForwardModeling(encsrc, dcal);
    matrix_transpose(&dcal[0], &trans[0], ng, nt);
    std::copy(pdata, pdata + numDataSamples, local_HOnA.getData() + i * numDataSamples);

    DEBUG() << format("parallel: sum HonA %.20f") % getSum(local_HOnA);

    TRACE() << "save the resd";
    const Matrix::value_type *obsDataBegin = local_D.getData() + i * numDataSamples;
    const Matrix::value_type *obsDataEnd   = obsDataBegin + numDataSamples;
    const Matrix::value_type *synDataBegin = local_HOnA.getData() + i * numDataSamples;
    resdSet[i] = variance(obsDataBegin, obsDataEnd, synDataBegin);

  }

  int count = local_D.getNumCol() * local_D.getNumRow();

  Matrix t4(N, N); /// HA' * U * SSqInv * U' * (D - HA)

	Matrix::value_type sum_local_D = pGetSum(local_D, N);
	Matrix::value_type sum_local_HOnA = pGetSum(local_HOnA, N);
  Matrix local_HA_Perturb(local_n, numDataSamples);
	pInitGamma(local_HOnA, local_HA_Perturb, nSamples);
	Matrix::value_type sum_HA_Pertrub = pGetSum(local_HA_Perturb, nSamples);
  Matrix local_perturbation(local_n, numDataSamples);
  pInitPerturbation(local_perturbation, local_HA_Perturb, rank, nSamples);
	Matrix::value_type sum_local_perturbation = pGetSum(local_perturbation, nSamples);
  std::transform(local_D.getData(), local_D.getData() + local_D.size(), local_perturbation.getData(), local_D.getData(), std::plus<Matrix::value_type>());
	Matrix::value_type sum_local_D2 = pGetSum(local_D, nSamples);
  Matrix local_gamma(local_n, numDataSamples);
  pInitGamma(local_perturbation, local_gamma, nSamples);
	Matrix::value_type sum_local_gamma = pGetSum(local_gamma, nSamples);
  Matrix local_band(local_n, numDataSamples);
  A_plus_B(local_HA_Perturb, local_gamma, local_band);
	Matrix::value_type sum_local_band = pGetSum(local_band, nSamples);
  Matrix local_matU(local_n, numDataSamples);
  Matrix local_matS(1, nSamples);
  Matrix local_matVt(local_n, nSamples);
	int info = pSvd(local_band, local_matU, local_matS, local_matVt, nSamples);
  Matrix matU2(N, numDataSamples);
  Matrix matVt2(N, nSamples);
	Matrix::value_type sum_local_matU = pGetSum(local_matU, nSamples);
  Matrix local_matSSqInv(N, N);
	if(	rank == 0) {
    DEBUG() << "parallel: sum of local_D: " << sum_local_D;
    DEBUG() << "parallel: sum of local_HOnA: " << sum_local_HOnA;
    DEBUG() << "parallel: sum of local_HA_Perturb: " << sum_HA_Pertrub;
    TRACE() << "parallel: initialize the perturbation";
    DEBUG() << "parallel: sum of local_perturbation: " << sum_local_perturbation;
    TRACE() << "parallel: add perturbation to observed data";
    DEBUG() << "parallel: sum of local_D with perturbation added: " << sum_local_D2;
    TRACE() << "parallel: calculate the gamma";
    DEBUG() << "parallel: sum of local_gamma: " << sum_local_gamma;
    TRACE() << "parallel: calculate the band";
    DEBUG() << "parallel: sum of local_band: " << sum_local_band;
    TRACE() << "parallel: svd of band";
    if ( info > 0 ) {
      ERROR() << "The algorithm computing SVD failed to converge.";
      exit( 1 );
    }
    DEBUG() << "parallel: sum of local_matU: " << sum_local_matU;
    DEBUG() << "parallel: sum of local_matS: " << getSum(local_matS);
    TRACE() << "parallel: calculate matSSqInv";
    int local_clip = clipPosition(local_matS);
    DEBUG() << "parallel: local_clip: " << local_clip;
    DEBUG() << "parallel: local_matS: ";
    local_matS.print();
    {
      Matrix::value_type *p = local_matSSqInv.getData();
      Matrix::value_type *s = local_matS.getData();
      const int nrow = local_matSSqInv.getNumRow();
      for (int i = 0; i < local_clip; i++) {
        p[i * nrow + i] = (1 / s[i]) * (1 / s[i]);
      }
    }
		//local_matSSqInv.print("local_matSSqInv.txt");
	}
  Matrix local_t0(local_n, numDataSamples); /// D - HA
	A_minus_B(local_D, local_HOnA, local_t0);
	Matrix::value_type sum_local_t0 = pGetSum(local_t0, nSamples);
  Matrix local_t1(local_n, nSamples); /// U' * (D - HA)
  pAlpha_ATrans_B_plus_beta_C(1.0, local_matU, 1, local_t0, 1, 0.0, local_t1, 1, nSamples);
	Matrix::value_type sum_local_t1 = pGetSum(local_t1, nSamples);
  Matrix local_t2(local_n, nSamples); /// U' * (D - HA)
  pAlpha_ATrans_B_plus_beta_C(1.0, local_HA_Perturb, 1, local_matU, 1, 0.0, local_t2, 1, nSamples);
	Matrix::value_type sum_local_t2 = pGetSum(local_t2, nSamples);
  Matrix local_t3(local_n, nSamples); /// U' * (D - HA)
  pAlpha_A_B_plus_beta_C(1.0, local_t2, 1, local_matSSqInv, 0, 0.0, local_t3, 1, nSamples);
	Matrix::value_type sum_local_t3 = pGetSum(local_t3, nSamples);
  Matrix local_t4(local_n, nSamples); /// U' * (D - HA)
  pAlpha_A_B_plus_beta_C(1.0, local_t3, 1, local_t1, 1, 0.0, local_t4, 1, nSamples);
	Matrix::value_type sum_local_t4 = pGetSum(local_t4, nSamples);
	/*
	char filename[20];
	sprintf(filename, "local_t4%d.txt", rank);
	local_t4.print(filename);
	*/
	if(rank == 0)
	{
    DEBUG() << "parallel: sum of matSSqInv: " << getSum(local_matSSqInv);
    TRACE() << "parallel: matSSQInv: ";
    DEBUG() << "parallel: sum of t0: " << sum_local_t0;
    DEBUG() << "parallel: sum of t1: " << sum_local_t1;
    DEBUG() << "parallel: sum of t2: " << sum_local_t2;
    DEBUG() << "parallel: sum of t3: " << sum_local_t3;
    DEBUG() << "parallel: sum of t4: " << sum_local_t4;
	}
  DEBUG() << "parallel: print HA' * U * SSqInv * U' * (D - HA)";
	local_t4.print();
	return local_t4;
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

  int local_n = velSet.size();
  std::vector<std::vector<int> > codes(local_n);

  int numDataSamples = ng * nt;
  std::vector<float> obsData(numDataSamples);
  std::vector<float> synData(numDataSamples);

  for (int i = 0; i < local_n; i++) {
    DEBUG() << format("init Lambda on velocity %2d/%d") % (i + 1) % velSet.size();

    ///  "making encoded shot, both sources and receivers";
    codes[i] = enkfRandomCodes.genPlus1Minus1(ns);

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
    assert(lambdaSet.getNumCol() == local_n);
    Matrix::value_type *plamda = lambdaSet.getData() + i * lambdaSet.getNumRow();
    const Matrix::value_type *pratio = ratioSet.getData() + i * ratioSet.getNumRow();
    plamda[0] = pratio[0] * resd / Wx2; /// lambdaX
    plamda[1] = pratio[1] * resd / Wz2; /// lamdaZ
    DEBUG() << "pratio[0] " << pratio[0];
    DEBUG() << "plambda[0] " << plamda[0];

  }

}

void EnkfAnalyze::pInitRatioPerturb(const Matrix& ratioSet, Matrix& ratioPerturb, int nsamples) const {
  pInitGamma(ratioSet, ratioPerturb, nsamples);
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

		//printf("%f ", avg);

    for (int icol = 0; icol < gamma.getNumCol(); icol++) {
      int idx = icol * gamma.getNumRow() + irow;
      gamma.getData()[idx] = p[idx] - avg;
    }
  }
	//printf("\n");
}

void EnkfAnalyze::pInitGamma(const Matrix& perturbation, Matrix& gamma, const int nSamples) const {
  assert(perturbation.isCompatible(gamma));

  const Matrix::value_type *p = perturbation.getData();
	std::vector<Matrix::value_type> sum_t(gamma.getNumRow(), 0.0);
	std::vector<Matrix::value_type> sum(gamma.getNumRow(), 0.0);

  for (int irow = 0; irow < gamma.getNumRow(); irow++) {
    for (int icol = 0; icol < gamma.getNumCol(); icol++) {
      sum_t[irow] += p[icol * gamma.getNumRow() + irow];
    }
	}

	MPI_Allreduce(&sum_t[0], &sum[0], gamma.getNumRow(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int irow = 0; irow < gamma.getNumRow(); irow++) {
		Matrix::value_type avg = sum[irow] / nSamples;
		//printf("%f ", avg);
    for (int icol = 0; icol < gamma.getNumCol(); icol++) {
      int idx = icol * gamma.getNumRow() + irow;
      gamma.getData()[idx] = p[idx] - avg;
    }
  }

}

void EnkfAnalyze::checkMatrix(const Matrix& a, const Matrix& b) const {
	for(int irow = 0 ; irow < a.getNumRow() ; irow ++)
	{
		for(int icol = 0 ; icol < a.getNumCol() ; icol ++)
		{
      int idx = icol * a.getNumRow() + irow;
			Matrix::value_type abserr = fabs(a.getData()[idx] - b.getData()[idx]);
			Matrix::value_type objerr = abserr / a.getData()[idx];
			if(objerr > 0.0001)
				printf("irow = %d, icol = %d, abserr = %lf, objerr = %lf\n", irow, icol, abserr, objerr);
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

std::vector<float> EnkfAnalyze::pCreateAMean(const std::vector<float*>& velSet, const int nSamples) const {
  int modelSize = fm.getnx() * fm.getnz();
	std::vector<float> ret(modelSize);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::vector<float> sum(modelSize);
  for (int i = 0; i < modelSize; i++) {
    for (size_t j = 0; j < velSet.size(); j++) {
      sum[i] += velSet[j][i];
    }
  }
	MPI_Allreduce(&sum[0], &ret[0], modelSize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i < modelSize; i++) {
		ret[i] /= nSamples;
	}

	return ret;
}

void EnkfAnalyze::check(std::vector<float> a, std::vector<float> b)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("Rank = %d, checking...\n", rank);
	int counter = 0;
	for(int i = 0 ; i < a.size() ; i ++)
	{
		float abserr = fabs(a[i] - b[i]);
		float relerr = abserr / a[i];
		if(relerr > 0.000001)
		{
			printf("Check index: %d %f %f %f\n", i, a[i], b[i], relerr);
			counter ++;
		}
		if(counter > 10)
			return;
	}
	if(counter == 0)
		printf("Pass!\n");
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

void EnkfAnalyze::pInitPerturbation(Matrix& perturbation, const Matrix &HA_Perturb, const int rank, const int nSamples) const {
  if (!initSigma) {
    initSigma = true;
		//how to make time seed in parallel program
    //int seed = 1;
    int seed = rank;
    double mean = 0;
    double maxHAP = std::abs(*std::max_element(HA_Perturb.getData(), HA_Perturb.getData() + HA_Perturb.size(), abs_less<float>));
		double totalMaxHAP = 0;
		MPI_Allreduce(&maxHAP, &totalMaxHAP, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		double sigma = initPerturbSigma(totalMaxHAP, sigmaFactor);
		DEBUG() << "parallel: sigmaIter0: " << sigma;
		sigmaIter0 = sigma;
		generator = new boost::variate_generator<boost::mt19937, boost::normal_distribution<> >(boost::mt19937(seed), boost::normal_distribution<>(mean, sigma));
  }

  std::generate(perturbation.getData(), perturbation.getData() + perturbation.size(), *generator);
}

void EnkfAnalyze::pInitPerturbation2(Matrix& perturbation, const Matrix &HA_Perturb, const int rank, const int nSamples) const {
  if (!initSigma) {
    initSigma = true;
    int seed = 1;
    double mean = 0;
    double maxHAP = std::abs(*std::max_element(HA_Perturb.getData(), HA_Perturb.getData() + HA_Perturb.size(), abs_less<float>));
		double totalMaxHAP = 0;
		MPI_Allreduce(&maxHAP, &totalMaxHAP, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		double sigma = initPerturbSigma(totalMaxHAP, sigmaFactor);
		DEBUG() << "parallel: sigmaIter0: " << sigma;
		sigmaIter0 = sigma;
		generator = new boost::variate_generator<boost::mt19937, boost::normal_distribution<> >(boost::mt19937(seed), boost::normal_distribution<>(mean, sigma));
  }

	std::vector<double> perturbation_t(nSamples * perturbation.getNumRow(), 0.0);
	if(rank == 0)
	{
		std::generate(&perturbation_t[0], &perturbation_t[0] + nSamples * perturbation.getNumRow(), *generator);
	}
	MPI_Scatter(&perturbation_t[0], perturbation.size(), MPI_DOUBLE, perturbation.getData(), perturbation.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
