#ifndef _YAFEL_SPARSE_COO_HPP
#define _YAFEL_SPARSE_COO_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_utils.hpp"
#include "lin_alg/sparse_matrix.hpp"
#include "lin_alg/construction_sparse_matrix.hpp"

#include <tuple>
#include <algorithm>
#include <cstdlib>
#include <cstdio>

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
  size_type nnz() {
    if(!_isCompressed) {
      compress();
    }
    return _nnz;
  }
  
  /*
   * Constructors
   */
  sparse_coo() : 
    _data(), _rows(0), _cols(0), 
    _nnz(0), _isCompressed(false)
  {}


  sparse_coo(const std::vector<triplet> &ts) : 
    _data(ts), _rows(0), _cols(0), 
    _nnz(0), _isCompressed(false) {
    
    compress();

    for(auto t : _data) {
      size_type i = std::get<0>(t);
      size_type j = std::get<1>(t);

      _rows = (_rows < i+1) ? i+1 : _cols;
      _cols = (_cols < j+1) ? j+1 : _cols;
    }
    
  }


  //sparse_coo(const sparse_coo &coo) : sparse_coo(coo._data) {}

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

      std::tuple<size_type, size_type> mid_ij(std::get<0>(_data[mid]), std::get<1>(_data[mid]));
      
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
  
  void add(size_type i, size_type j, value_type val) {
      _data.push_back(triplet(i,j,val));

    // update some matrix properties
    _isCompressed = false;
    _rows = (_rows < i+1) ? i+1 : _rows;
    _cols = (_cols < j+1) ? j+1 : _cols;
  }

  void preallocate(size_type N) {
    _data.reserve(N);
  }

  // utility routine, useful for testing if nothing else
  bool operator==(sparse_coo<dataType> &rhs) {
    compress();
    rhs.compress();
    
    if(_rows != rhs.rows() ||
       _cols != rhs.cols() ||
       _nnz != rhs.nnz()) {
      
      return false;
    }

    for(size_type i=0; i<_data.size(); ++i) {
      if(_data[i] != rhs._data[i]) {
	return false;
      }
    }

    return true;
  }

  // Return a constant reference to the triplets, to be used in construction 
  // of sparsity patterns and other sparse matrix formats.
  // There is no need to make full copies for this, since the sparsity 
  // pattern needs to allocate its own memory anyway.
  const std::vector<triplet> & get_triplets() {
    compress();
    return _data;
  }
  
private:
  std::vector<triplet> _data;
  size_type _rows;
  size_type _cols;
  size_type _nnz;
  bool _isCompressed;

  void compress() {
    
    if(_isCompressed) {
      return;
    }
    
    // Sort data using standard tuple operator< : sorts first by i, then j, then val.
    std::sort(_data.begin(), _data.end());

    
    // Compress data by accumulating all entries with identical (i,j) into a single entry
    size_type current_location = 0;
    for(size_type i=1; i<_data.size(); ++i) {
      if(std::get<0>(_data[current_location])==std::get<0>(_data[i]) &&
	 std::get<1>(_data[current_location])==std::get<1>(_data[i])) {

	std::get<2>(_data[current_location]) += std::get<2>(_data[i]);
      }
      else {
	_data[++current_location] = _data[i];
      }
    }

    _nnz = current_location+1;
    _isCompressed = true;
    _data.resize(_nnz);
  }
};

YAFEL_NAMESPACE_CLOSE

#endif
