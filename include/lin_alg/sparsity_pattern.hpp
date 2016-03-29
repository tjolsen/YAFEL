#ifndef _YAFEL_SPARSITY_PATTERN_HPP
#define _YAFEL_SPARSITY_PATTERN_HPP

/*
 * Base class to define interface for the sparsity pattern for sparse matrices.
 * This will be specialized to implement CSR, CSC, BCSR, and BCSC sparse matrices
 */

#include <vector>


YAFEL_NAMESPACE_OPEN

template<typename T>
class sparsity_pattern {

public:
  typedef std::vector<std::size_t> container_type;
  typedef typename container_type::size_type size_type;
  
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<const T&>(*this);}
};


YAFEL_NAMESPACE_CLOSE

#endif
