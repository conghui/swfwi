/*
 * Matrix.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rice
 */

#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <cstring>
#include <mpi.h>

#include "Matrix.h"
#include "logger.h"


namespace {

std::streampos getFileSize(const std::string &filePath) {
  std::streampos fsize = 0;
  std::ifstream file(filePath.c_str(), std::ios::binary);

  if (!file) {
    ERROR() << "cannot open file: " << filePath;
    exit(EXIT_FAILURE);
  }

  fsize = file.tellg();
  file.seekg( 0, std::ios::end );
  fsize = file.tellg() - fsize;
  file.close();

  return fsize;
}

}

Matrix::Matrix(int ncol, int nrow) :
  mData(NULL), mNumRow(nrow), mNumCol(ncol) {
  mData = (value_type *)malloc(nrow * ncol * sizeof(value_type));
  if (mData == NULL) {
    ERROR() << __PRETTY_FUNCTION__ << ": malloc fails";
    perror("reason is:");
    exit(0);
  }
  std::fill(mData, mData + nrow * ncol, 0);
}

Matrix::~Matrix() {
  free(mData);
}


void Matrix::readFromFile(const std::string &filename) {
  std::ifstream ifs(filename.c_str());
  assert(ifs.good());

  size_t fileSize = getFileSize(filename);
  assert(fileSize == mNumRow * mNumCol * sizeof(value_type));

  ifs.read(reinterpret_cast<char *>(&mData[0]), fileSize);

  ifs.close();
}

void Matrix::print() const {
  std::stringstream ss;
  TRACE() << "print in column-major order";
  for (int col = 0; col < getNumCol(); col++) {
    for (int row = 0; row < getNumRow(); row++) {
      ss << std::setprecision(4) << std::fixed << mData[col * mNumRow + row] << " ";
    }
    ss << "\n\n";
  }

  DEBUG() << ss.str();
}

//row-majored print
void Matrix::print(char *filename) const {
	FILE *f = fopen(filename, "w");
  for (int row = 0; row < getNumRow(); row++) {
		for (int col = 0; col < getNumCol(); col++) {
			fprintf(f, "%lf ", mData[col * mNumRow + row]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void Matrix::printInfo(char *filename) const {
	FILE *f = fopen(filename, "w");
	fprintf(f, "%d %d", getNumRow(), getNumCol());
	fclose(f);
}

bool Matrix::isCompatible(const Matrix &rhs) const {
  return mNumRow == rhs.mNumRow && mNumCol == rhs.mNumCol;
}

int Matrix::getNumCol() const {
  return mNumCol;
}

int Matrix::getNumRow() const {
  return mNumRow;
}

int Matrix::size() const {
  return getNumRow() * getNumCol();
}


const Matrix::value_type *Matrix::getData() const {
  return &mData[0];
}

Matrix::value_type *Matrix::getData() {
  return const_cast<value_type *>((static_cast<const Matrix *>(this))->getData());
}


void A_plus_B(const Matrix &A, const Matrix &B, Matrix &C) {
  assert(A.isCompatible(B));
  assert(A.isCompatible(C));
  std::transform(A.getData(), A.getData() + A.size(), B.getData(), C.getData(), std::plus<double>());
}

Matrix::value_type getSum(const Matrix &M) {
  const Matrix::value_type *p = M.getData();
  const int size = M.getNumRow() * M.getNumCol();
  Matrix::value_type sum = std::accumulate(p, p + size, 0.0);
  return sum;
}

Matrix::value_type pGetSum(const Matrix &M, const int nSamples) {
  const Matrix::value_type *p = M.getData();
  const int size = M.getNumRow() * M.getNumCol();
  Matrix::value_type sum = std::accumulate(p, p + size, 0.0);
	//printf("sum = %e\n", sum);
	Matrix::value_type ret = 0;
	MPI_Reduce(&sum, &ret, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//printf("ret = %e\n", ret);
  return ret;
}

Matrix::value_type pGetSum2(const Matrix &M, const int nSamples) {
	Matrix total_M(nSamples, M.getNumRow());
  int size = M.getNumRow() * M.getNumCol();
	MPI_Gather(M.getData(), size, MPI_DOUBLE, total_M.getData(), size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  const Matrix::value_type *p = total_M.getData();
	size = total_M.getNumRow() * total_M.getNumCol();
  Matrix::value_type sum = std::accumulate(p, p + size, 0.0);
  return sum;
}

void A_minus_B(const Matrix &A, const Matrix &B, Matrix &C) {
  assert(A.isCompatible(B));
  assert(A.isCompatible(C));
  std::transform(A.getData(), A.getData() + A.size(), B.getData(), C.getData(), std::minus<Matrix::value_type>());
}

int clipPosition(const Matrix &M) {
  const Matrix::value_type relativeError = 0.001f;
  const int N = M.getNumRow();
  int pos;
  Matrix::value_type sum = 0;
  const Matrix::value_type *m = M.getData();

  for (pos = 0; pos < N - 1; pos++) {
    sum += std::abs(m[pos]);
    if (relativeError / (1 - relativeError) * sum > (N - pos + 1) * m[pos + 1]) {
      break;
    }
  }

  return pos + 1;
}

float getSum(std::vector<float> &v)
{
	float sum = 0.0;
	for(int i = 0 ; i < v.size() ; i ++)
		sum += v[i];
	return sum;
}
