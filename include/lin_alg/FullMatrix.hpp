#ifndef _YAFEL_FULLMATRIX_HPP
#define _YAFEL_FULLMATRIX_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

class FullMatrix {
  
private:
  unsigned rows;
  unsigned cols;
  bool transposed;
  double *data;
  inline unsigned indexOf(unsigned i, unsigned j) const;
  FullMatrix inverse2x2() const;
  FullMatrix inverse3x3() const;
  
public:
  FullMatrix(unsigned m, unsigned n);
  FullMatrix(unsigned m, unsigned n, double val);
  FullMatrix(unsigned n);
  FullMatrix(unsigned n, double val);
  FullMatrix(const FullMatrix & src);
  FullMatrix(const sparse_csr & csr);
  FullMatrix(const sparse_coo & coo);
  ~FullMatrix();
  
  double & operator()(unsigned i, unsigned j) const;
  FullMatrix operator*(const FullMatrix & rhs) const;
  FullMatrix &operator+=(const FullMatrix &rhs);
  FullMatrix operator+(const FullMatrix &rhs) const;
  FullMatrix &operator-=(const FullMatrix &rhs);
  FullMatrix operator-(const FullMatrix &rhs) const;
  Vector operator*(const Vector & rhs) const;
  inline unsigned getRows() const { return (!transposed) ? rows : cols; }
  inline unsigned getCols() const { return (!transposed) ? cols : rows; }
  double det() const;
  FullMatrix getTransposed() const;
  FullMatrix getInverse() const;
  Vector slice_col(unsigned col) const;
  Vector slice_row(unsigned row) const;
  void transpose();
  void print();
};

YAFEL_NAMESPACE_CLOSE

#endif
