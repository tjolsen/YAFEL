#ifndef __YAFEL_SPARSE_CSR_HPP
#define __YAFEL_SPARSE_CSR_HPP

#include "lin_alg/sparse_matrix.hpp"
#include "lin_alg/access_sparse_matrix.hpp"
#include "lin_alg/csr_sparsity_pattern.hpp"

#include <iostream>
#include <algorithm>
#include <cstdlib>

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class sparse_csr : public access_sparse_matrix<sparse_csr<dataType>, dataType> {

public:
  using container_type =  typename sparse_matrix<sparse_csr<dataType>, dataType>::container_type;
  using size_type =  typename sparse_matrix<sparse_csr<dataType>, dataType>::size_type;
  using value_type = typename sparse_matrix<sparse_csr<dataType>, dataType>::value_type;
  using reference =  typename sparse_matrix<sparse_csr<dataType>, dataType>::reference;
  using triplet = typename sparse_matrix<sparse_csr<dataType>, dataType>::triplet;

  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}
  size_type nnz() const {return _data.size();}

  sparse_csr(const std::vector<triplet> &ts) : _zero(0) {

    
    // handle empty matrix. Currently, it'll break as soon as operator() is called.
    if(ts.size() == 0) {
      return;
    }


    // std::less sorts triplets in (row, col, val) form
    // into an order that is useful for construction of
    // sparse_csr matrices.
    std::vector<triplet> triplets(ts);
    std::sort(triplets.begin(), triplets.end());
    
    //get the number of rows+1
    row_ptr.resize(std::get<0>(triplets[triplets.size()-1])+2);

    size_type curr_row = std::get<0>(triplets[0]);
    size_type curr_col = std::get<1>(triplets[0]);

    //initialize matrix size. will grow when new triplets are processed
    _rows = curr_row+1;
    _cols = curr_col+1;

    for(size_type r=0; r<=curr_row; ++r) {
      row_ptr[r] = 0;
    }
    col_index.push_back(curr_col);
    _data.push_back(std::get<2>(triplets[0]));
    size_type idx = 0;
    
    for(size_type i=1; i<triplets.size(); ++i) {
      size_type row = std::get<0>(triplets[i]);
      size_type col = std::get<1>(triplets[i]);
      dataType val = std::get<2>(triplets[i]);

      _rows = ((row+1) > _rows) ? row+1 : _rows;
      _cols = ((col+1) > _cols) ? col+1 : _cols;
      
      if(row == curr_row && col==curr_col) {
	_data[idx] += val;
      }
      else {
	++idx;
	if(row != curr_row) {
	  for(size_type r=curr_row+1; r<=row; ++r) {
	    row_ptr[r] = idx;
	  }
	  curr_row = row;
	}
	col_index.push_back(col);
	curr_col = col;
	_data.push_back(val);
      }
    }

    row_ptr[_rows] = _data.size();
  }

  template<typename T>
  sparse_csr(sparse_matrix<T,dataType> &sp) : sparse_csr(sp.get_triplets()) {}


  /*
   * Index operators
   * Defining _OPTIMIZED turns off bounds checking
   */
  value_type operator()(size_type i, size_type j) const {
#ifndef _OPTIMIZED
    if(i>=_rows || j >= cols) {
      std::cerr << "error:sparse_csr:operator() out of bounds\n";
      exit(1);
    }
#endif
    bool in_sparsity = true;
    size_type idx = index_of(i,j,in_sparsity);
    return (in_sparsity) ? _data[idx] : dataType(0);
  }

  value_type& operator()(size_type i, size_type j) {
#ifndef _OPTIMIZED
    if(i>=_rows || j >= cols) {
      std::cerr << "error:sparse_csr:operator() out of bounds\n";
      exit(1);
    }
#endif

    bool in_sparsity = true;
    size_type idx = index_of(i,j,in_sparsity);
    return (in_sparsity) ? _data[idx] : _zero;
  }

  /*
   * Functions to get sparsity patterns.
   * csr_sparsity_pattern_copy contains a deep copy of row_ptr, col_index
   * csr_sparsity_pattern_reference contains const references to these vectors
   */
  csr_sparsity_pattern_copy get_sparsity_pattern_copy() const {
    return csr_sparsity_pattern_copy(row_ptr, col_index);
  }

  csr_sparsity_pattern_reference get_sparsity_pattern_reference() const {
    return csr_sparsity_pattern_reference(row_ptr, col_index);
  }


  // easier to write algorithms if these are public
  std::vector<size_type> row_ptr;
  std::vector<size_type> col_index;
  container_type _data;

private:
  size_type _rows;
  size_type _cols;
  value_type _zero; // <-- this holds a zero so that I can return value_type& from a operator() call
  
  inline size_type index_of(size_type i, size_type j, bool &in_sparsity) {
    for(size_type idx=row_ptr[i]; idx<row_ptr[i+1]; ++idx) {
      if(j==col_index[idx]) {
	return idx;
      }
    }
    
    in_sparsity=false;
    return _data.size(); // this can't possibly point to anything, so easy to trigger errors
  }

}; //end class sparse_csr

YAFEL_NAMESPACE_CLOSE

#endif
