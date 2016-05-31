/*
 * pMatrix.h
 *
 *  Created on: May 19, 2016
 *      Author: Bingwei Chen
 */
#ifndef H_P_MATRIX
#define H_P_MATRIX

#include<stdlib.h>
#include<memory.h>
#include<stdio.h>
#include<math.h>
#include<iostream>
#include "Matrix.h"
#include "mpi.h"
extern "C" {
int blacs_get_(int *val, int *what, int *icontxt);
int	blacs_gridinit_(int *icontxt, char *order, int *nprow, int *npcol);
int	blacs_gridinfo_(int *icontxt, int *nprow, int *npcol, int *myrow, int *mycol);
int	blacs_barrier_(int *icontxt, char *scope);
int	blacs_exit_(int *i_continue);
int numroc_(int *m, int *mb, int *myrow, int *i_zero, int *nprow );
int	descinit_(int *descA, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *icontxt, int *lld, int *info);
void pdgeadd_(char *trans, int *m, int *n, double *alpha, double *A, int *IA, int *JA, int *descA, double *beta, double *C, int *IC, int *JC, int *descC);
int	pdgesvd_(char *JOBU, char *JOBV, int *m, int *n, double *A, int *IA, int *JA, int *descA, double *S, double *U, int *IU, int *JU, int *descU, double *Vt, int *IVt, int *JVt, int *descVt, double *work, int *lwork, int *info);
void pdgemm_(char* TRANSA, char* TRANSB, int * M, int * N, int * K, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB, double * BETA, double * C, int * IC, int * JC, int * DESCC);
}

class pMatrix
{
	public:
		pMatrix(double *_data, int _grow, int _gcol, int _lrow, int _lcol, bool _global);
		pMatrix(int _grow, int _gcol, int _lrow, int _lcol, bool _global);
		~pMatrix();
		void initGrid();
		void initContent();
		void print();
		void print(const char filename[]);
		void printInfo(const char head[]);
		void read(char *filename);
		void test();
		double* getData();
		int* getDesc();
		int getGRow();
		int getGCol();
		int getMp();
		int getNq();
		void setMp(int _mp);
		void setNq(int _nq);
		static void init(int proSize);
		static void finalize();
		static int ictxt;	//content

	private:
		static int nprow;	//row processes number in grid
		static int npcol;	//column processes number in grid
		static int myrow;	//this process row index in grid
		static int mycol;	//this process column index in grid
		int mb;	//matrix row partitioning
		int nb;	//matrix col partitioning
		double *data;	//data array, column major!!!!!!
		int grow;	//global row of matrix
		int gcol;	//global col of matrix
		int lrow;	//real row of matrix
		int lcol;	//real col of matrix
		int mp;	//row needed in every process
		int nq;	//column needed in every process
		bool global;	//whether the matrix is a global matrix or just a part
		int *desc;	//matrix description array
		int lld;	//leading dimension of matrix
		int info;		//information
};

class pMatrixMM
{
	public:
		//A: M*K, B: K*N, C: M*N
		pMatrixMM(char _transa, char _transb, int _M, int _N, int _K, double _alpha, pMatrix *_A, pMatrix *_B, double _beta, pMatrix *_C);
		void run();
		void row_col(int &M, int &N, int &grow, int &gcol, int &lrow, int &lcol, bool &global, int size, int kind);
	private:
		char transa;
		char transb;
		pMatrix *A;
		pMatrix *B;
		pMatrix *C;
		int M;
		int N;
		int K;
		double alpha;
		double beta;
};

class pMatrixSVD
{
	public:
		pMatrixSVD(pMatrix *_A, pMatrix *_U, pMatrix *_S, pMatrix *_Vt);
		~pMatrixSVD();
		void run();
		int getInfo();
	private:
		char JOBU;
		char JOBV;
		pMatrix *A;
		pMatrix *U;
		pMatrix *S;
		pMatrix *Vt;
		int M;
		int N;
		double *work;
		int lwork;
		int info;
};

void pAlpha_A_B_plus_beta_C(double alpha, Matrix &A, int kindA, Matrix &B, int kindB, double beta, Matrix &C, int kindC, const int nSamples);
void pAlpha_ATrans_B_plus_beta_C(double alpha, Matrix &tA, int kindA, Matrix &tB, int kindB, double beta, Matrix &tC, int kindC, const int nSamples);
void pAlpha_A_B_plus_beta_C(char transa, char transb, double alpha, Matrix &A, int kindA, Matrix &B, int kindB, double beta, Matrix &C, int kindC, const int nSamples);
int pSvd(Matrix &A, Matrix &U, Matrix &S, Matrix &Vt, const int nSamples);

#endif
