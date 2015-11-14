#ifndef __YAFEL_SPARSE_CSR_HPP
#define __YAFEL_SPARSE_CSR_HPP

#include "lin_alg/sparse_matrix.hpp"
#include "lin_alg/csr_sparsity_pattern.hpp"

#include <algorithm>

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class sparse_csr : public sparse_matrix<sparse_csr<dataType>, dataType> {

public:
  using container_type =  typename sparse_matrix<sparse_csr<dataType>, dataType>::container_type;
  using size_type =  typename sparse_matrix<sparse_csr<dataType>, dataType>::size_type;
  using value_type = typename sparse_matrix<sparse_csr<dataType>, dataType>::value_type;
  using reference =  typename sparse_matrix<sparse_csr<dataType>, dataType>::reference;
  using triplet = typename sparse_matrix<sparse_csr<dataType>, dataType>::triplet;

  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}
  size_type nnz() const {return _data.size();}

  sparse_csr(std::vector<triplet> triplets) {
    // std::less sorts triplets in (row, col, val) form
    // into an order that is useful for construction of
    // sparse_csr matrices.
    std::sort<triplet,std::less>(triplets.begin(), triplets.end());
    
    //get the number of rows+1
    row_ptr.resize(std::get<0>(*triplets.end())+1);


    size_type curr_row = std::get<0>(triplets[0]);
    size_type curr_col = std::get<1>(triplets[0]);
    dataType curr_val =  std::get<2>(triplets[0]);

    //initialize matrix size. will grow when new triplets are processed
    _rows = curr_row;
    _cols = curr_col;

    row_ptr[0] = 0;
    col_index.push_back(curr_col);
    
    for(size_type i=1; i<triplets.size(); ++i) {
      size_type row = std::get<0>(triplets[i]);
      size_type col = std::get<1>(triplets[i]);
      dataType val = std::get<2>(triplets[i]);

      _rows = ((row+1) > _rows) ? row+1 : _rows;
      _cols = ((col+1) > _cols) ? col+1 : _cols;

      if(row == curr_row && col==curr_col) {
	// compress new value into accumulator
	curr_val += val;
      }
      else if(row == curr_row) {
	// store accumulated value, move to next col
	_data.push_back(curr_val);
	curr_val = val;

	col_index.push_back(curr_col);
	curr_col = col;
      }
    }
    
  }

  template<typename T>
  sparse_csr(sparse_matrix<T,dataType> &sp) : sparse_csr(sp.get_triplets()) {}


  value_type operator()(size_type i, size_type j) const {
    bool in_sparsity = false;
    size_type idx = index_of(i,j,in_sparsity);
    return (in_sparsity) ? _data[idx] : _zero;
  }

  value_type& operator()(size_type i, size_type j) {
    bool in_sparsity = false;
    size_type idx = index_of(i,j,in_sparsity);
    return (in_sparsity) ? _data[idx] : _zero;
  }

private:
  std::vector<size_type> row_ptr;
  std::vector<size_type> col_index;
  container_type _data;
  size_type _rows;
  size_type _cols;
  value_type _zero; // <-- this holds a zero so that I can return value_type& from a operator() call

  inline size_type index_of(size_type i, size_type j, bool &in_sparsity) const {
    in_sparsity = true;
    size_type lo, hi, mid;
    do {
      lo = row_ptr[i];
      hi = row_ptr[i+1];
      mid = lo + (hi-lo)/2;
      
      size_type jmid = col_index[mid];
      if( jmid < j) {
	lo = mid+1;
      }
      else if(jmid > j) {
	hi = mid;
      }
      else {
	in_sparsity = true;
	break;
      }
    } while (lo < hi);
    
    return mid;
  }

}; //end class sparse_csr

YAFEL_NAMESPACE_CLOSE

#endif
