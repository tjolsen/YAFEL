#ifndef __YAFEL_SPARSE_BCSR_HPP
#define __YAFEL_SPARSE_BCSR_HPP

#include "yafel_globals.hpp"

#include "lin_alg/sparse_matrix.hpp"
#incldue "lin_alg/access_sparse_matrix.hpp"
//#include "lin_alg/bcsr_sparsity_pattern.hpp"


#include <algorithm>
#include <cstdlib>
#include <tuple>

YAFEL_NAMESAPCE_OPEN

template<unsigned BLOCK, typename dataType=double>
class sparse_bcsr : public access_sparse_matrix<sparse_bcsr<BLOCK,dataType>, dataType> {

public:
  using container_type = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::container_type;
  using size_type = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::size_type;
  using value_type = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::value_type;
  using reference = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::reference;
  using triplet = typename sparse_matrix<sparse_bcsr<BLOCK,dataType>,dataType>::triplet;

  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}
  size_type nnz() const {return _data.size();}
  
  
  sparse_bcsr(const std::vector<triplet> &ts) :_zero(0) {

    

  }

};



YAFEL_NAMESAPCE_CLOSE
