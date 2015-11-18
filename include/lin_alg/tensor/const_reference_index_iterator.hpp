#ifndef __YAFEL_CONST_REFERENCE_INDEX_ITERATOR_HPP
#define __YAFEL_CONST_REFERENCE_INDEX_ITERATOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/generic_index_iterator.hpp"

YAFEL_NAMESPACE_OPEN

template<typename TENSTYPE, unsigned DIM, unsigned RANK, int... IDXS>
class const_reference_index_iterator : public generic_index_iterator<DIM,RANK,IDXS...> {

private:
  const TENSTYPE &_t;

  template<int ...S>
  typename TENSTYPE::value_type apply_indices(const seq<S...> &) const {
    return _t(std::get<S>(this->indices)...);
  }
  
public:
  template<typename ...Args>
  const_reference_index_iterator(const TENSTYPE &t, Args... _indices) : 
    generic_index_iterator<DIM,RANK,IDXS...>(_indices...), 
    _t(t)
  {
    static_assert(sizeof...(_indices)==RANK,
		  "const_reference_index_iterator: wrong numbero f indices");
  }
  
  typename TENSTYPE::value_type operator*() const {
    return apply_indices(this->sequence);
  }
};


YAFEL_NAMESPACE_CLOSE

#endif
