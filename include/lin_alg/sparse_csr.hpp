#ifndef __YAFEL_SPARSE_CSR_HPP
#define __YAFEL_SPARSE_CSR_HPP

#include "lin_alg/sparse_matrix.hpp"
#include "lin_alg/csr_sparsity_pattern.hpp"


template<typename dataType>
class sparse_csr : public sparse_matrix<sparse_csr<dataType>, dataType> {

public:
  typedef typename sparse_matrix<sparse_csr<dataType>, dataType>::container_type container_type;
  typedef typename sparse_matrix<sparse_csr<dataType>, dataType>::size_type size_type;
  typedef typename sparse_matrix<sparse_csr<dataType>, dataType>::value_type value_type;
  typedef typename sparse_matrix<sparse_csr<dataType>, dataType>::reference reference;
  typedef typename sparse_matrix<sparse_csr<dataType>, dataType>::triplet triplet;

  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}
  size_type nnz() const {return data.size();}


  sparse_csr(std::vector<triplet> &triplets) {
    // std::less sorts triplets in (row, col, val) form
    // into an order that is useful for construction of
    // sparse_csr matrices.
    std::sort(triplets.begin(), triplets.end(), std::less);
    
    //get the number of rows+1
    row_ptr.resize(std::get<0>(*triplets.end())+1);


    size_type curr_row = std::get<0>(triplets[0]);
    size_type curr_col = std::get<1>(triplets[0]);
    dataType curr_val =  std::get<2>(triplets[0]);

    row_ptr[0] = 0;
    col_index.push_back(curr_col);
    
    for(size_type i=1; i<triplets.size(); ++i) {
      size_type row = std::get<0>(triplets[i]);
      size_type col = std::get<1>(triplets[i]);
      dataType val = std::get<2>(triplets[i]);
      
      if(row == curr_row && col==curr_col) {
	// compress new value into accumulator
	curr_val += val;
      }
      else if(row == curr_row) {
	// store accumulated value, move to next col
	
	curr_col = col;
	
      }
    }
    
  }



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
  container_type data;
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

};


#endif
