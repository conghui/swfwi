/*
 * pMatrix.cpp
 *
 *  Created on: May 19, 2016
 *      Author: Bingwei Chen
 */
#include "pMatrix.h"

int i_zero = 0;
int i_one = 1;
int i_negone = -1;
double d_zero = 0;
double d_one = 1;
double d_negone = -1;
int i_continue = 0;
char c_scope = 'A';
int DLEN_ = 9;

int pMatrix::nprow;	//row processes number in grid
int pMatrix::npcol;	//column processes number in grid
int pMatrix::myrow;	//this process row index in grid
int pMatrix::mycol;	//this process column index in grid
int pMatrix::ictxt;	//content

void pMatrix::test()
{
	printf("myrow = %d, mycol = %d\n", myrow, mycol);
}

void pMatrix::init(int proSize)
{
	nprow = 1;
	npcol = proSize;
	//get content
	blacs_get_( &i_negone, &i_zero, &ictxt );
	//initial grid
	blacs_gridinit_( &ictxt, "C", &nprow, &npcol );
	//get rank
	blacs_gridinfo_( &ictxt, &nprow, &npcol, &myrow, &mycol);
}

void pMatrix::finalize()
{
	blacs_exit_(&i_zero);
}

int pMatrix::getMp()
{
	return mp;
}

int pMatrix::getNq()
{
	return nq;
}

void pMatrix::setMp(int _mp)
{
	mp = _mp;
	descinit_(desc, &grow, &gcol, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
}

void pMatrix::setNq(int _nq)
{
	nq = _nq;
	descinit_(desc, &grow, &gcol, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
}

//If matrix is the whole matrix, grow == lrow, gcol == lcol. If the matrix is one part of the whole matrix, grow and gcol is the row and column number of the whole matrix, lrow and lcol is the real row and column number
pMatrix::pMatrix(int _grow, int _gcol, int _lrow, int _lcol, bool _global)
{
	grow = _grow;
	gcol = _gcol;
	lrow = _lrow;
	lcol = _lcol;
	global = _global;
	desc = (int *)malloc(sizeof(int) * DLEN_);
	if(global)
	{
		if(mycol == 0)
			data = (double *)malloc(sizeof(double) * grow * gcol);
		else
			data = (double *)malloc(sizeof(double));
	}
	else
		data = (double *)malloc(sizeof(double) * lrow * lcol);
	initGrid();
}

pMatrix::pMatrix(double *_data, int _grow, int _gcol, int _lrow, int _lcol, bool _global)
{
	data = _data;
	grow = _grow;
	gcol = _gcol;
	lrow = _lrow;
	lcol = _lcol;
	global = _global;
	desc = (int *)malloc(sizeof(int) * DLEN_);
	int rank;
	initGrid();
}

void pMatrix::initGrid()
{
	//mb = lrow;
	nb = lcol;
	mb = nb;	//in svd, mb should be equal to nb
	//printf("nprow = %d, npcol = %d, mb = %d, nb = %d\n", nprow, npcol, mb, nb);
	mp = numroc_(&grow, &mb, &myrow, &i_zero, &nprow);
	nq = numroc_(&gcol, &nb, &mycol, &i_zero, &npcol);
	lld = std::max(mp, 1);
	descinit_(desc, &grow, &gcol, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
}

void pMatrix::printInfo(const char head[])
{
	printf("%s: myrow = %d, mycol = %d, grow = %d, gcol = %d, lrow = %d, lcol = %d, mb = %d, nb = %d, mp = %d, nq = %d\n", head, myrow, mycol, grow, gcol, lrow, lcol, mb, nb, mp, nq);
}

void pMatrix::print()
{
	//printf("myrow = %d, mycol = %d\n", myrow, mycol);
	for(int i = 0 ; i < mp ; i ++)
	{
		for(int j = 0 ; j < nq ; j ++)
			printf("%lf ", data[j * mp + i]);
		printf("\n");
	}
	printf("\n");
}

void pMatrix::print(const char filename[])
{
	FILE *f = fopen(filename, "w");
	fprintf(f, "myrow = %d, mycol = %d\n", myrow, mycol);
	for(int i = 0 ; i < mp ; i ++)
	{
		for(int j = 0 ; j < nq ; j ++)
			fprintf(f, "%lf ", data[j * mp + i]);
		fprintf(f, "\n");
	}
	printf("\n");
	fclose(f);
}

void pMatrix::read(char *filename)
{
	if(global && mycol != 0)
		return;
	//printf("row = %d, col = %d\n", lrow, lcol);
	FILE *f = fopen(filename,"r");
	for(int i = 0 ; i < mp ; i ++)
	{
		for(int j = 0 ; j < nq ; j ++)
			fscanf(f, "%lf", &data[j * mp + i]);
	}
	fclose(f);
}

double* pMatrix::getData()
{
	return data;
}

int* pMatrix::getDesc()
{
	return desc;
}

int pMatrix::getGRow()
{
	return grow;
}

int pMatrix::getGCol()
{
	return gcol;
}

pMatrix::~pMatrix()
{
	free(desc);
}

pMatrixMM::pMatrixMM(char _transa, char _transb, int _M, int _N, int _K, double _alpha, pMatrix *_A, pMatrix *_B, double _beta, pMatrix *_C)
{
	transa = _transa;
	transb = _transb;
	A = _A;
	B = _B;
	C = _C;
	M = _M;
	N = _N;
	K = _K;
	alpha = _alpha;
	beta = _beta;
}

void pMatrixMM::run()
{
	pdgemm_(&transa, &transb, &M, &N, &K, &alpha, A->getData(), &i_one, &i_one, A->getDesc(), B->getData(), &i_one, &i_one, B->getDesc(), &beta, C->getData(), &i_one, &i_one, C->getDesc());
}

pMatrixSVD::pMatrixSVD(pMatrix *_A, pMatrix *_U, pMatrix *_S, pMatrix *_Vt)
{
	JOBU = 'V';
	JOBV = 'V';
	A = _A;
	U = _U;
	S = _S;
	Vt = _Vt;
	M = A->getGRow();
	N = A->getGCol();
	lwork = -1;
	work = (double *)malloc(sizeof(double));
}

pMatrixSVD::~pMatrixSVD()
{
	free(work);
}

void pMatrixSVD::run()
{
	pdgesvd_(&JOBU, &JOBV, &M, &N, A->getData(), &i_one, &i_one, A->getDesc(), S->getData(), U->getData(), &i_one, &i_one, U->getDesc(), Vt->getData(), &i_one, &i_one, Vt->getDesc(), work, &lwork, &info);
	lwork = work[0];
	free(work);
	work = (double *)malloc(sizeof(double) * lwork);
	pdgesvd_(&JOBU, &JOBV, &M, &N, A->getData(), &i_one, &i_one, A->getDesc(), S->getData(), U->getData(), &i_one, &i_one, U->getDesc(), Vt->getData(), &i_one, &i_one, Vt->getDesc(), work, &lwork, &info);
}

int pMatrixSVD::getInfo()
{
	return info;
}

void pAlpha_A_B_plus_beta_C(double alpha, Matrix &tA, int kindA, Matrix &tB, int kindB, double beta, Matrix &tC, int kindC, const int nSamples)
{
	pAlpha_A_B_plus_beta_C('N', 'N', alpha, tA, kindA, tB, kindB, beta, tC, kindC, nSamples);
}

void pAlpha_ATrans_B_plus_beta_C(double alpha, Matrix &tA, int kindA, Matrix &tB, int kindB, double beta, Matrix &tC, int kindC, const int nSamples)
{
	pAlpha_A_B_plus_beta_C('T', 'N', alpha, tA, kindA, tB, kindB, beta, tC, kindC, nSamples);
}

void row_col(int &M, int &N, int &grow, int &gcol, int &lrow, int &lcol, bool &global, int size, int kind)
{
	if(kind == 0)
	{
		MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	if(kind / 2 == 1)
		M *= size;
	if(kind % 2 == 1)
		N *= size;
	grow = M;
	gcol = N;
	lrow = M;
	lcol = N;
	global = true;

	if(kind / 2 == 1)
	{
		lrow /= size;
		global = false;
	}
	if(kind % 2 == 1)
	{
		lcol /= size;
		global = false;
	}
}

//kind: 0 stands for global, 1 stands for column partitioning, 2 stands for row partitioning
void pAlpha_A_B_plus_beta_C(char transa, char transb, double alpha, Matrix &tA, int kindA, Matrix &tB, int kindB, double beta, Matrix &tC, int kindC, const int nSamples)
{
	int rank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int grow_A, gcol_A, lrow_A, lcol_A;
	bool global_A;
	int M_A = tA.getNumRow();
	int N_A = tA.getNumCol();
	row_col(M_A, N_A, grow_A, gcol_A, lrow_A, lcol_A, global_A, size, kindA);

	int grow_B, gcol_B, lrow_B, lcol_B;
	bool global_B;
	int M_B = tB.getNumRow();
	int N_B = tB.getNumCol();
	row_col(M_B, N_B, grow_B, gcol_B, lrow_B, lcol_B, global_B, size, kindB);

	int grow_C, gcol_C, lrow_C, lcol_C;
	bool global_C;
	int M_C = tC.getNumRow();
	int N_C = tC.getNumCol();
	row_col(M_C, N_C, grow_C, gcol_C, lrow_C, lcol_C, global_C, size, kindC);

	int M, N, K;
	if(transa == 'N')
	{
		M = grow_A;
		K = gcol_A;
	}
	else
	{
		M = gcol_A;
		K = grow_A;
	}
	if(transb == 'N')
		N = gcol_B;
	else
		N = grow_B;

	pMatrix::init(size);
	pMatrix A(tA.getData(), grow_A, gcol_A, lrow_A, lcol_A, global_A);
	pMatrix B(tB.getData(), grow_B, gcol_B, lrow_B, lcol_B, global_B);
	pMatrix C(tC.getData(), grow_C, gcol_C, lrow_C, lcol_C, global_C);

	pMatrixMM mm(transa, transb, M, N, K, alpha, &A, &B, beta, &C);
	mm.run();
}

int pSvd(Matrix &tA, Matrix &tU, Matrix &tS, Matrix &tVt, const int nSamples)
{
	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int M = tA.getNumRow();
	int N = nSamples;
	int minSize = std::min(M, N);

	int grow_band = M;
	int gcol_band = N;
	int lrow_band = M;
	int lcol_band = N / size;

	int grow_U = M;
	int gcol_U = minSize;
	int lrow_U = M;
	int lcol_U = minSize / size;

	int grow_S = 1;
	int gcol_S = minSize;
	int lrow_S = 1;
	int lcol_S = minSize;

	int grow_Vt = minSize;
	int gcol_Vt = N;
	int lrow_Vt = minSize;
	int lcol_Vt = N / size;

	pMatrix::init(size);
	pMatrix band(tA.getData(), grow_band, gcol_band, lrow_band, lcol_band, false);
	pMatrix U(tU.getData(), grow_U, gcol_U, lrow_U, lcol_U, false);
	pMatrix S(tS.getData(), grow_S, gcol_S, lrow_S, lcol_S, false);
	pMatrix Vt(tVt.getData(), grow_Vt, gcol_Vt, lrow_Vt, lcol_Vt, false);

	/*
	char filename[20];
	sprintf(filename, "output_band%d", rank);
	band.print(filename);
	*/

	pMatrixSVD svd(&band, &U, &S, &Vt);
	svd.run();
	return svd.getInfo();

	/*
	sprintf(filename, "U%d.txt", rank);
	U.print(filename);
	sprintf(filename, "S%d.txt", rank);
	S.print(filename);
	sprintf(filename, "Vt%d.txt", rank);
	Vt.print(filename);
	*/
}
