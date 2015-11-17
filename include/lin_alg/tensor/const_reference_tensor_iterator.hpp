#ifndef _CONST_REFERENCE_TENSOR_ITERATOR_HPP
#define _CONST_REFERENCE_TENSOR_ITERATOR_HPP

#include "generic_tensor_iterator.hpp"

template<typename TENSTYPE, unsigned DIM, unsigned RANK>
class const_reference_tensor_iterator : public generic_tensor_iterator<DIM,RANK> {

private:
  const TENSTYPE &_t;

  template<int ...S>
  typename TENSTYPE::value_type apply_indices(const seq<S...> &) const{
    return _t(std::get<S>(this->indices)...);
  }
  
public:
  const_reference_tensor_iterator(const TENSTYPE &t) : _t(t)
  {}
  
  typename TENSTYPE::value_type operator*() const {
    return apply_indices(this->sequence);
  }
};

#endif
