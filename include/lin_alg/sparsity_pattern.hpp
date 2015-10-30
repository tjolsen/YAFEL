#ifndef __YAFEL_SPARSITY_PATTERN_HPP
#define __YAFEL_SPARSITY_PATTERN_HPP

/*
 * Base class to define interface for the sparsity pattern for sparse matrices.
 * This will be specialized to implement CSR, CSC, BCSR, and BCSC sparse matrices
 */

#include <vector>

template<typename T>
class sparsity_pattern {

public:
  typedef std::vector<std::size_t> container_type;
  typedef typename container_type::size_type size_type;
  
  size_type operator()(size_type i, size_type j, bool &in_sparsity) const {
    return static_cast<const T&>(*this)(i,j,in_sparsity);
  }
  
  size_type rows() const {return static_cast<const T&>(*this).rows();}
  size_type cols() const {return static_cast<const T&>(*this).cols();}
  size_type nnz() const {return static_cast<const T&>(*this).nnz();}

  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<const T&>(*this);}
};

#endif
