/*
 * Matrix.h
 *
 *  Created on: Dec 4, 2014
 *      Author: rice
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <string>

class Matrix {
 public:
  typedef double value_type;
  Matrix(int ncol, int nrow);
  ~Matrix();
  void readFromFile(const std::string &filename);
  void print() const;

  bool isCompatible(const Matrix &rhs) const;
  int getNumCol() const;
  int getNumRow() const;
  int size() const;
  value_type *getData();
  const value_type *getData() const;

 private:
  value_type *mData;
  int mNumRow;
  int mNumCol;
};

//void extract_each_col_from_each_row_mean(const Matrix &A, Matrix &ret);

//void alpha_A_B_plus_beta_C(const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, double alpha, const Matrix &A, const Matrix &B, double beta, Matrix &C);
void alpha_ATrans_B_plus_beta_C(double alpha, const Matrix &A, const Matrix &B, double beta, Matrix &C);
void alpha_A_B_plus_beta_C(double alpha, const Matrix &A, const Matrix &B, double beta, Matrix &C);
void A_plus_B(const Matrix &A, const Matrix &B, Matrix &C);
void A_minus_B(const Matrix &A, const Matrix &B, Matrix &C);

Matrix::value_type getSum(const Matrix &M);

/// compute clip position
int clipPosition(const Matrix &M);

#endif /* MATRIX_H_ */
