#ifndef _YAFEL_SPARSE_COO_HPP
#define _YAFEL_SPARSE_COO_HPP

#include <stdio.h>
#include <stdlib.h>
#include "yafel_globals.hpp"
#include "lin_alg/sparse_csr.hpp"

#define INSERTION_SORT_THRESHOLD 10

YAFEL_NAMESPACE_OPEN

class sparse_coo {
  
  friend class sparse_csr;
  
private:
  static const unsigned default_capacity = 10;
  unsigned rows;
  unsigned cols;
  unsigned size;
  unsigned capacity;
  unsigned oldNNZ;
  bool isSorted;
  bool consistentNNZ;
  unsigned *row_index;
  unsigned *col_index;
  double *data;
  void init();
  void mergeSort(unsigned start, unsigned end, unsigned *Brow, unsigned *Bcol, double *Bval);
  void merge(unsigned start, unsigned middle, unsigned end, unsigned *Brow, unsigned *Bcol, double *Bval);
  void insertionSort(unsigned start, unsigned end);
  
public:
  sparse_coo();
  sparse_coo(unsigned N);
  sparse_coo(unsigned m, unsigned n);
  sparse_coo(const sparse_coo & src);
  sparse_coo(const sparse_csr & csr);
  ~sparse_coo();
  sparse_coo & operator=(const sparse_coo &rhs);
  inline unsigned getRows() const {return rows;}
  inline unsigned getCols() const {return cols;}
  inline unsigned getSize() const {return size;}
  inline bool isSquare() const {return rows==cols;}
  unsigned nnz();
  void add(unsigned row, unsigned col, double val);
  void sort();
  void compress();
  void print_sparse();
  void load_from_file(const char *fname);
};

YAFEL_NAMESPACE_CLOSE

#endif
