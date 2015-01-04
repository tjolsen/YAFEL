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
  static const int default_capacity = 10;
  int rows;
  int cols;
  int size;
  int capacity;
  int oldNNZ;
  bool isSorted;
  bool consistentNNZ;
  int *row_index;
  int *col_index;
  double *data;
  void init();
  void mergeSort(int start, int end, int *Brow, int *Bcol, double *Bval);
  void merge(int start, int middle, int end, int *Brow, int *Bcol, double *Bval);
  void insertionSort(int start, int end);
  
public:
  sparse_coo();
  sparse_coo(int N);
  sparse_coo(int m, int n);
  sparse_coo(const sparse_coo & src);
  sparse_coo(const sparse_csr & csr);
  ~sparse_coo();
  sparse_coo & operator=(const sparse_coo &rhs);
  int getRows() const {return rows;}
  int getCols() const {return cols;}
  int getSize() const {return size;}
  bool isSquare() const {return rows==cols;}
  int nnz();
  void add(int row, int col, double val);
  void sort();
  void compress();
  void print_sparse();
  void load_from_file(const char *fname);
};

YAFEL_NAMESPACE_CLOSE

#endif
