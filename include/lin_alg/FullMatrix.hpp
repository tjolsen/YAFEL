#ifndef _YAFEL_FULLMATRIX_HPP
#define _YAFEL_FULLMATRIX_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

class FullMatrix {
  
private:
  int rows;
  int cols;
  bool transposed;
  double *data;
  inline int indexOf(int i, int j) const;
  FullMatrix inverse2x2() const;
  FullMatrix inverse3x3() const;
  
public:
  FullMatrix(int m, int n);
  FullMatrix(int m, int n, double val);
  FullMatrix(int n);
  FullMatrix(int n, double val);
  FullMatrix(const FullMatrix & src);
  FullMatrix(const sparse_csr & csr);
  FullMatrix(const sparse_coo & coo);
  ~FullMatrix();
  
  double & operator()(int i, int j) const;
  FullMatrix operator*(const FullMatrix & rhs) const;
  Vector operator*(const Vector & rhs) const;
  int getRows() const { return rows; }
  int getCols() const { return cols; }
  double det() const;
  FullMatrix getTransposed() const;
  FullMatrix getInverse() const;
  Vector slice_col(int col) const;
  Vector slice_row(int row) const;
  void transpose();
  void print();
};

YAFEL_NAMESPACE_CLOSE

#endif
