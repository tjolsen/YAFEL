#ifndef _YAFEL_REFERENCE_TENSOR_ITERATOR_HPP
#define _YAFEL_REFERENCE_TENSOR_ITERATOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/generic_tensor_iterator.hpp"

YAFEL_NAMESPACE_OPEN

template<typename TENSTYPE, unsigned DIM, unsigned RANK>
class reference_tensor_iterator : public generic_tensor_iterator<DIM,RANK> {

private:
  TENSTYPE &_t;

  template<int ...S>
  typename TENSTYPE::value_type & apply_indices(seq<S...> &) {
    return _t(std::get<S>(this->indices)...);
  }

public:
  reference_tensor_iterator(TENSTYPE &t) : _t(t)
  {}

  typename TENSTYPE::value_type & operator*() {
    return apply_indices(this->sequence);
  }
};

YAFEL_NAMESPACE_CLOSE

#endif
