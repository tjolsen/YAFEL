#ifndef _YAFEL_SPARSE_COO_HPP
#define _YAFEL_SPARSE_COO_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_matrix.hpp"
#include "lin_alg/construction_sparse_matrix.hpp"

#include <tuple>
#include <cstdlib>
#include <algorithm>

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class sparse_coo : public construction_sparse_matrix<sparse_coo<dataType>, dataType> {
public:
  typedef typename construction_sparse_matrix<sparse_coo<dataType>,dataType>::size_type size_type;
  typedef typename construction_sparse_matrix<sparse_coo<dataType>,dataType>::value_type value_type;
  typedef typename construction_sparse_matrix<sparse_coo<dataType>,dataType>::reference reference;
  typedef typename construction_sparse_matrix<sparse_coo<dataType>,dataType>::triplet triplet;


  size_type rows() const {return _rows;}
  size_type cols() const {return _cols;}

  /*
   * Will implement 2 kinds of operator() functions with different const-ness.
   * The non-const one will be faster for repeated use (which you shouldn't really be doing, btw),
   * as it will be permitted to call the compress() function and find the desired element via
   * binary search.
   *
   * The const version will not be permitted to do this, so it must do a linear scan through ALL
   * triples to ensure proper aggregation of the value.
   */
  value_type operator()(size_type i, size_type j) const {
    value_type val(0);
    
    for( const triplet& t : _data) {
      if(std::get<0>(t)==i && std::get<1>(t)==j) {
	val += std::get<2>(t);
      }
    }
    
    return val;
  }

  value_type operator()(size_type i, size_type j) {
    if(!_isCompressed) {
      compress();
    }
    std::tuple<size_type, size_type> target(i,j);
    size_type lo(0), hi(_nnz), mid(0);
    do {
      if(hi==lo) {
	mid = lo;
	break;
      }
      mid = lo + (hi-lo)/2;

      std::tuple<size_type, size_type> mid_ij(std::get<0>(_data[mid]), std:;get<1>(_data[mid]));
      
      if(mid_ij < target) {
	hi=mid;
      }
      else if(mid_ij > target) {
	lo = mid+1;
      }
      else {
	break;
      }
    } while(true);
    
    triplet& t = _data[mid];
    if(std::get<0>(t)==i && std::get<1>(t)==j) {
      return std::get<2>(t);
    }
    else {
      return 0;
    }
  }

private:
  std::vector<triplet> _data;
  size_type _rows;
  size_type _cols;
  size_type _nnz;
  bool _isCompressed;

  void compress();
};

YAFEL_NAMESPACE_CLOSE

#endif
