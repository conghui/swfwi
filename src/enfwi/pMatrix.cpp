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
	initGrid();
}

void pMatrix::initGrid()
{
	//mb = lrow;
	nb = lcol;
	mb = nb;	//in svd, mb should be equal to nb
	printf("nprow = %d, npcol = %d, mb = %d, nb = %d\n", nprow, npcol, mb, nb);
	mp = numroc_(&grow, &mb, &myrow, &i_zero, &nprow);
	nq = numroc_(&gcol, &nb, &mycol, &i_zero, &npcol);
	lld = std::max(mp, 1);
	descinit_(desc, &grow, &gcol, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info);
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
	printf("a\n");
	if(global && mycol != 0)
		return;
	printf("row = %d, col = %d\n", lrow, lcol);
	FILE *f = fopen(filename,"r");
	for(int i = 0 ; i < mp ; i ++)
	{
		for(int j = 0 ; j < nq ; j ++)
			fscanf(f, "%lf", &data[j * mp + i]);
	}
	fclose(f);
	printf("b\n");
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
}

pMatrixMM::pMatrixMM(double _alpha, pMatrix *_A, pMatrix *_B, double _beta, pMatrix *_C)
{
	transa = 'N';
	transb = 'N';
	A = _A;
	B = _B;
	C = _C;
	M = A->getGRow();
	N = B->getGCol();
	K = A->getGCol();
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

void pMatrixSVD::run()
{
	pdgesvd_(&JOBU, &JOBV, &M, &N, A->getData(), &i_one, &i_one, A->getDesc(), S->getData(), U->getData(), &i_one, &i_one, U->getDesc(), Vt->getData(), &i_one, &i_one, Vt->getDesc(), work, &lwork, &info);
	lwork = work[0];
	free(work);
	work = (double *)malloc(sizeof(double) * lwork);
	pdgesvd_(&JOBU, &JOBV, &M, &N, A->getData(), &i_one, &i_one, A->getDesc(), S->getData(), U->getData(), &i_one, &i_one, U->getDesc(), Vt->getData(), &i_one, &i_one, Vt->getDesc(), work, &lwork, &info);
	//printf("info = %d\n", info);
}

void pAlpha_A_B_plus_beta_C(double alpha, Matrix &tA, Matrix &tB, double beta, Matrix &tC, const int nSamples)
{
	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int M = tA.getNumRow();
	int K = nSamples;
	int N = nSamples;

	int grow_A = M;
	int gcol_A = K;
	int lrow_A = M;
	int lcol_A = K / size;

	int grow_B = K;
	int gcol_B = N;
	int lrow_B = K;
	int lcol_B = N;

	int grow_C = M;
	int gcol_C = N;
	int lrow_C = M;
	int lcol_C = N / size;

	pMatrix::init(size);
	pMatrix A(tA.getData(), grow_A, gcol_A, lrow_A, lcol_A, false);
	pMatrix B(tB.getData(), grow_B, gcol_B, lrow_B, lcol_B, true);
	pMatrix C(tC.getData(), grow_C, gcol_C, lrow_C, lcol_C, false);

	pMatrixMM mm(alpha, &A, &B, beta, &C);
	mm.run();

	char scope = 'A';
	blacs_barrier_(&pMatrix::ictxt, &scope);
	/*
	char filename[20];
	sprintf(filename, "local_Perturb%d.txt", rank);
	A.print(filename);
	sprintf(filename, "local_Matrix%d.txt", rank);
	B.print(filename);
	sprintf(filename, "local_t5%d.txt", rank);
	C.print(filename);
	*/
}

void pSvd(Matrix &tA, Matrix &tU, Matrix &tS, Matrix &tVt, const int nSamples)
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

	char filename[20];
	sprintf(filename, "output_band%d", rank);
	band.print(filename);

	printf("here1\n");
	pMatrixSVD svd(&band, &U, &S, &Vt);
	svd.run();

	printf("here2\n");

	sprintf(filename, "U%d.txt", rank);
	U.print(filename);
	sprintf(filename, "S%d.txt", rank);
	S.print(filename);
	sprintf(filename, "Vt%d.txt", rank);
	Vt.print(filename);
	MPI_Barrier(MPI_COMM_WORLD);

}
