#ifndef _YAFEL_SPARSE_CSR_HPP
#define _YAFEL_SPARSE_CSR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

class sparse_coo; //forward declaration. cyclic dependancy
class DirBC;

class sparse_csr {
  
  friend class sparse_coo;
  friend class DirBC;

private:
  unsigned rows;
  unsigned cols;
  unsigned nnz;
  unsigned *row_ptr;
  unsigned *col_index;
  double *data;
  
  
public:
  sparse_csr(sparse_coo & src);
  sparse_csr(const sparse_csr & src);
  ~sparse_csr();
  void init_from_coo(sparse_coo & src);
  inline unsigned getRows() const {return rows;}
  inline unsigned getCols() const {return cols;}
  inline unsigned getNNZ() const { return nnz;}
  double operator()(unsigned i, unsigned j, bool & flag) const;
  double operator()(unsigned i, unsigned j) const;
  sparse_csr & operator*=(double a);
  sparse_csr operator*(const sparse_csr & rhs) const;
  Vector operator*(const Vector & rhs) const;
  Vector slice_col(unsigned col) const;
  Vector slice_row(unsigned row) const;
  void zero_col(unsigned col);
  void zero_row(unsigned row);
  void assign(unsigned row, unsigned col, double val);
  void print_sparse();
};

YAFEL_NAMESPACE_CLOSE

#endif
