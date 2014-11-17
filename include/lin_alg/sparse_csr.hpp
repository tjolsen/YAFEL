#ifndef _YAFEL_SPARSE_CSR_HPP
#define _YAFEL_SPARSE_CSR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

class sparse_coo; //forward declaration. cyclic dependancy

class sparse_csr {
  
  friend class sparse_coo;
  
private:
  int rows;
  int cols;
  int nnz;
  int *row_ptr;
  int *col_index;
  double *data;
  
  
public:
  sparse_csr(sparse_coo & src);
  sparse_csr(const sparse_csr & src);
  ~sparse_csr();
  void init_from_coo(sparse_coo & src);
  int getRows() const {return rows;}
  int getCols() const {return cols;}
  int getNNZ() const { return nnz;}
  double operator()(int i, int j, bool & flag) const;
  double operator()(int i, int j) const;
  sparse_csr operator*(const sparse_csr & rhs) const;
  Vector operator*(const Vector & rhs) const;
  Vector slice_col(int col) const;
  Vector slice_row(int row) const;
  void zero_col(int col);
  void zero_row(int row);
  void assign(int row, int col, double val);
  void print_sparse();
};

YAFEL_NAMESPACE_CLOSE

#endif
