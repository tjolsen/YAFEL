#ifndef __YAFEL_SPARSE_BCSR_HPP
#define __YAFEL_SPARSE_BCSR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_matrix.hpp"
#include "lin_alg/access_sparse_matrix.hpp"
//#include "lin_alg/bcsr_sparsity_pattern.hpp"

#include <iostream> // <-- for debugging only
#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <cassert>

YAFEL_NAMESPACE_OPEN

template<unsigned BLOCK, typename dataType=double>
class sparse_bcsr : public access_sparse_matrix<sparse_bcsr<BLOCK,dataType>, dataType> {

public:
  using container_type = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::container_type;
  using size_type = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::size_type;
  using value_type = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::value_type;
  using reference = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::reference;
  using triplet = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::triplet;

  size_type rows() const {return _brows*BLOCK;}
  size_type cols() const {return _bcols*BLOCK;}
  size_type nnz() const {return data.size();}
  
  std::vector<size_type> brow_ptr;
  std::vector<size_type> bcol_index;
  container_type data;

  struct bcsr_compare {
    bool operator()(const triplet &t1, const triplet &t2) {
      size_type I1, I2, J1, J2, i1, i2, j1, j2;
      I1 = std::get<0>(t1)/BLOCK;
      I2 = std::get<0>(t2)/BLOCK;
      J1 = std::get<1>(t1)/BLOCK;
      J2 = std::get<1>(t2)/BLOCK;
      i1 = std::get<0>(t1)%BLOCK;
      i2 = std::get<0>(t2)%BLOCK;
      j1 = std::get<1>(t1)%BLOCK;
      j2 = std::get<1>(t2)%BLOCK;
      
      if(I1 < I2) {
	return true;
      }
      else {
	if(I1 > I2) return false;
	
	if(J1 < J2) {
	  return true;
	}
	else {
	  if(J1>J2) return false;

	  if(i1<i2) {
	    return true;
	  }
	  else {
	    if (i1>i2) return false;

	    if(j1<j2) {
	      return true;
	    }
	    else {
	      return false;
	    }
	  }
	}
      }
    }
  };


  /*
   * Construction from other sparse_matrix objects (including, perhaps especially, sparse_coo)
   *
   */
  template<typename T>
  sparse_bcsr(sparse_matrix<T,dataType> &sp) : sparse_bcsr(sp.get_triplets()) {}
  

  /*
   * Construction from triplets.
   * This is the main workhorse of the data structure, taking an unordered (potentially uncompressed)
   * vector of triplets and performing the required sorting to convert it into a sparse_bcsr ordering.
   * 
   * This function DOES MODIFY the input triplet vector "ts" (by sorting)
   */
  sparse_bcsr(const std::vector<triplet> & _ts) : _zero(0) {
    std::vector<triplet> ts(_ts);
    block_stride = BLOCK*BLOCK;
    assert(ts.size() > 0 && "Doesn't support empty matrices (yet)");

    //sort the triplets into bcsr order using the bcsr_compare struct
    std::sort(ts.begin(), ts.end(), bcsr_compare());
    
    size_type brows = std::get<0>(*ts.rbegin())/BLOCK+1;
    brow_ptr.resize(brows+1, size_type(0));

    size_type curr_brow = std::get<0>(ts[0])/BLOCK;
    size_type curr_bcol = std::get<1>(ts[0])/BLOCK;
  
    // count number of blocks (single pass thru tuples. saves mem alloc time from push_back calls)
    size_type nblocks = 1;
    for(size_type i=1; i<ts.size(); ++i) {
      if(curr_brow == std::get<0>(ts[i])/BLOCK &&
	 curr_bcol == std::get<1>(ts[i])/BLOCK) {
      }
      else {
	++nblocks;
      }
      curr_brow = std::get<0>(ts[i])/BLOCK;
      curr_bcol = std::get<1>(ts[i])/BLOCK;
    }
    
    // resize storage containers
    bcol_index.resize(nblocks, size_type(0));
    data.resize(nblocks*block_stride, value_type(0));

    //construct from triplets
    curr_brow = std::get<0>(ts[0])/BLOCK;
    curr_bcol = std::get<1>(ts[0])/BLOCK;
    
    //_rows = std::get<0>(ts[0])+1;
    //_cols = std::get<1>(ts[0])+1;
    
    for(size_type r=0; r<=curr_brow; ++r) {
      brow_ptr[r] = 0;
    }
    
    bcol_index[0] = curr_bcol;
    data[index(0, std::get<0>(ts[0])%BLOCK, std::get<1>(ts[0])%BLOCK)] = std::get<2>(ts[0]);
    
    size_type bidx = 0;
    for(size_type idx=1; idx<ts.size(); ++idx) {
      
      size_type I = std::get<0>(ts[idx])/BLOCK;
      size_type J = std::get<1>(ts[idx])/BLOCK;
      size_type i = std::get<0>(ts[idx])%BLOCK;
      size_type j = std::get<1>(ts[idx])%BLOCK;
      value_type val = std::get<2>(ts[idx]);
      

      //keep track of number of brows, bcols. 
      // not yet sure if I should enforce matrix sizes to be multiples of BLOCK
      _brows = (I+1)>_brows ? I+1 : _brows;
      _bcols = (J+1)>_bcols ? J+1 : _bcols;
      
      
      if(I==curr_brow && J==curr_bcol) {
	data[index(bidx, i, j)] += val;
      }
      else {
	++bidx;

	if(curr_brow != I) {
	  for(size_type r=curr_brow+1; r<=I; ++r) {
	    brow_ptr[r] = bidx;
	  }
	  
	  curr_brow = I;
	}
	curr_bcol = J;

	bcol_index[bidx] = J;
	data[index(bidx, i, j)] += val;
      }
    }
    
    brow_ptr[_brows] = nblocks;
  }


  value_type operator()(size_type i, size_type j) const {
#ifndef _OPTIMIZED
    if(i>= _brows*BLOCK || j >=_bcols*BLOCK) {
      std::cerr << "error:sparse_bcsr:operator() out of bounds\n";
      exit(1);
    }
#endif
    
    bool in_sparsity = true;
    size_type bidx = bindex_of(i/BLOCK, j/BLOCK, in_sparsity);
    
    return (in_sparsity) ? data[index(bidx, i%BLOCK, j%BLOCK)] : dataType(0);
  }

private:
  size_type _brows;
  size_type _bcols;
  size_type _nnz;
  size_type block_stride;
  value_type _zero;

  inline size_type index(size_type bidx, size_type i, size_type j) const {
    return bidx*block_stride + i*BLOCK + j;
  }

  inline size_type bindex_of(size_type I, size_type J, bool &in_sparsity) const {
    for(size_type idx=brow_ptr[I]; idx<brow_ptr[I+1]; ++idx) {
      if(J == bcol_index[idx]) {
        in_sparsity = true;
        return idx;
      }
    }
    in_sparsity=false;
    return data.size();
  }


}; // end class sparse_bcsr



YAFEL_NAMESPACE_CLOSE

#endif
