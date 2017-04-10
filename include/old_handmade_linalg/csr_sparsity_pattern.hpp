#ifndef _YAFEL_CSR_SPARSITY_PATTERN_HPP
#define _YAFEL_CSR_SPARSITY_PATTERN_HPP


#include "sparsity_pattern.hpp"
#include "sparse_matrix.hpp"

#include <algorithm>

YAFEL_NAMESPACE_OPEN

class csr_sparsity_pattern_copy : public sparsity_pattern<csr_sparsity_pattern_copy> {

public:
  using container_type = typename sparsity_pattern<csr_sparsity_pattern_copy>::container_type;
  using size_type = typename sparsity_pattern<csr_sparsity_pattern_copy>::size_type;
  
  container_type row_ptr;
  container_type col_index;

  csr_sparsity_pattern_copy(const container_type &rp, const container_type &ci) :
    row_ptr(rp), col_index(ci)
  {}
  
  csr_sparsity_pattern_copy(container_type &&rp, container_type &&ci) :
    row_ptr(rp), col_index(ci)
  {}

};


class csr_sparsity_pattern_reference : public sparsity_pattern<csr_sparsity_pattern_reference> {

public:
  using container_type = typename sparsity_pattern<csr_sparsity_pattern_reference>::container_type;
  using size_type = typename sparsity_pattern<csr_sparsity_pattern_reference>::size_type;
  
  const container_type & row_ptr;
  const container_type & col_index;

  csr_sparsity_pattern_reference(const container_type &rp, const container_type &ci) :
    row_ptr(rp), col_index(ci)
  {}
  
  csr_sparsity_pattern_reference(container_type &&rp, container_type &&ci) :
    row_ptr(rp), col_index(ci)
  {}

};



YAFEL_NAMESPACE_CLOSE

#endif
