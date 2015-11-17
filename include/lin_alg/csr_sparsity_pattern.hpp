#ifndef __YAFEL_CSR_SPARSITY_PATTERN_HPP
#define __YAFEL_CSR_SPARSITY_PATTERN_HPP


#include "lin_alg/sparsity_pattern.hpp"
#include "lin_alg/sparse_matrix.hpp"

#include <algorithm>

YAFEL_NAMESPACE_OPEN

class csr_sparsity_pattern_copy : public sparsity_pattern<csr_sparsity_pattern> {

public:
  typedef typename sparsity_pattern<csr_sparsity_pattern>::container_type container_type;
  typedef typename sparsity_pattern<csr_sparsity_pattern>::size_type size_type;
  
  container_type row_ptr;
  container_type col_index;


  csr_sparsity_pattern(const container_type &rp, const container_type &ci) :
    row_ptr(rp), col_index(ci)
  {}
  
  csr_sparsity_pattern(container_type &&rp, container_type &&ci) :
    row_ptr(rp), col_index(ci)
  {}

};


class csr_sparsity_pattern_reference : public sparsity_pattern<csr_sparsity_pattern> {

public:
  typedef typename sparsity_pattern<csr_sparsity_pattern>::container_type container_type;
  typedef typename sparsity_pattern<csr_sparsity_pattern>::size_type size_type;
  
  const container_type & row_ptr;
  const container_type & col_index;

  csr_sparsity_pattern(const container_type &rp, const container_type &ci) :
    row_ptr(rp), col_index(ci)
  {}
  
  csr_sparsity_pattern(container_type &&rp, container_type &&ci) :
    row_ptr(rp), col_index(ci)
  {}

};



YAFEL_NAMESPACE_CLOSE

#endif
