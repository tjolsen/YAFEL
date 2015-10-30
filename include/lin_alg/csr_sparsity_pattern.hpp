#ifndef __YAFEL_CSR_SPARSITY_PATTERN_HPP
#define __YAFEL_CSR_SPARSITY_PATTERN_HPP


#include "lin_alg/sparsity_pattern.hpp"

YAFEL_NAMESPACE_OPEN

class csr_sparsity_pattern : public sparsity_pattern<csr_sparsity_pattern> {

public:
  typedef typename sparsity_pattern<csr_sparsity_pattern>::container_type container_type;
  typedef typename sparsity_pattern<csr_sparsity_pattern>::size_type size_type;

  container_type row_ptr;
  container_type col_index;

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
