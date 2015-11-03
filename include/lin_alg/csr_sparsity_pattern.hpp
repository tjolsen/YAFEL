#ifndef __YAFEL_CSR_SPARSITY_PATTERN_HPP
#define __YAFEL_CSR_SPARSITY_PATTERN_HPP


#include "lin_alg/sparsity_pattern.hpp"
#include "lin_alg/sparse_matrix.hpp"

#include <algorithm>

YAFEL_NAMESPACE_OPEN

class csr_sparsity_pattern : public sparsity_pattern<csr_sparsity_pattern> {

public:
  typedef typename sparsity_pattern<csr_sparsity_pattern>::container_type container_type;
  typedef typename sparsity_pattern<csr_sparsity_pattern>::size_type size_type;
  
  container_type row_ptr;
  container_type col_index;


  /*
   * Constructors
   */
  
  // copy const references
  csr_sparsity_pattern(const container_type &rp, const container_type &ci) :
    row_ptr(rp), col_index(ci)
  {}
  
  // move ctor for row_ptr and col_index
  csr_sparsity_pattern(container_type &&rp, container_type &&ci) :
    row_ptr(rp), col_index(ci)
  {}

  
  // construct from a vector of triplets
  template<typename dataType>
  csr_sparsity_pattern(std::vector< std::tuple<size_type, size_type, dataType> > & triplets) {
    std::sort(triplets.begin(), triplets.end(), std::less);
    
    
    
  }
  
  
  size_type operator()(size_type i, size_type j, bool &in_sparsity) {
    in_sparsity = false;

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
	lo = mid;
	in_sparsity = true;
	break;
      }
    } while (lo < hi);
    
    return mid;
  }
};

YAFEL_NAMESPACE_CLOSE

#endif
