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
  void print(char *filename) const;
  void printInfo(char *filename) const;

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

void alpha_A_B_plus_beta_C(double alpha, Matrix &A, Matrix &B, double beta, Matrix &C);
void alpha_ATrans_B_plus_beta_C(double alpha, Matrix& A, Matrix& B, double beta, Matrix& C);
void A_plus_B(const Matrix &A, const Matrix &B, Matrix &C);
void A_minus_B(const Matrix &A, const Matrix &B, Matrix &C);

Matrix::value_type getSum(const Matrix &M);
Matrix::value_type pGetSum(const Matrix &M, const int nSamples);
Matrix::value_type pGetSum2(const Matrix &M, const int nSamples);

/// compute clip position
int clipPosition(const Matrix &M);

#endif /* MATRIX_H_ */
