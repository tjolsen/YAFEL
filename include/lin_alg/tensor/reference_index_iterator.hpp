#ifndef _YAFEL_REFERENCE_INDEX_ITERATOR_HPP
#define _YAFEL_REFERENCE_INDEX_ITERATOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/generic_index_iterator.hpp"

YAFEL_NAMESPACE_OPEN

template<typename TENSTYPE, unsigned DIM, unsigned RANK, int... IDXS>
class reference_index_iterator : public generic_index_iterator<DIM,RANK,IDXS...> {

private:
  TENSTYPE &_t;

  template<int ...S>
  typename TENSTYPE::value_type & apply_indices(seq<S...> &) {
    return _t(std::get<S>(this->indices)...);
  }
  
public:
  template<typename ...Args>
  reference_index_iterator(TENSTYPE &t, Args... _indices) : 
    generic_index_iterator<DIM,RANK,IDXS...>(_indices...), 
    _t(t)
  {}
  
  typename TENSTYPE::value_type & operator*() {
    return apply_indices(this->sequence);
  }
};

YAFEL_NAMESPACE_CLOSE

#endif
