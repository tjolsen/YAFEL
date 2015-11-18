#ifndef __YAFEL_TENSOR_HPP
#define __YAFEL_TENSOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"


YAFEL_NAMESPACE_OPEN

template<unsigned DIM, unsigned RANK, typename dataType=double>
class Tensor : public TensorExpression< Tensor<DIM,RANK,dataType>, DIM, RANK, dataType> {
public:
  using value_type = typename TensorExpression<Tensor<DIM,RANK,dataType>,DIM,RANK,dataType>::value_type;
  
  template<typename ...Args>
  value_type & operator()(Args ...args) {return _data[index(args...)];}

  template<typename ...Args>
  value_type operator()(Args ...args) const {return _data[index(args...)];}
  
  //get a reference_tensor_iterator to the beginning
  reference_tensor_iterator<Tensor<DIM,RANK,dataType>,DIM,RANK> begin() {
    return reference_tensor_iterator<Tensor<DIM,RANK,dataType>,DIM,RANK>(*this);
  }
  
  //initialize a zero-tensor
  Tensor() {
    for(unsigned i=0; i<_tensorStorage(DIM,RANK); ++i) {
      _data[i] = 0;
    }
  }

  // copy a tensor or tensor-expression into a new (explicitly stored) tensor
  template <typename T>
  Tensor(const TensorExpression<T,DIM,RANK,dataType> &tens) {

    auto tens_it = tens.begin();

    for(auto TI = reference_tensor_iterator<Tensor<DIM,RANK,dataType>,DIM,RANK>(*this); 
	!TI.end(); TI.next(), tens_it.next()) {
      *TI = *tens_it;
    }
  }

private:
  value_type _data[_tensorStorage(DIM,RANK)];

  template<typename T>
  inline T index(T i) const {
    return i;
  }

  template<typename T, typename ...Args>
  inline T index(T i, Args... args) const {
    return DIM*index(args...) + i;
  }

  template<typename T, int ...S, typename ...Args>
  void copy_element(const seq<S...> &, 
		    const std::tuple<Args...> &indices, 
		    const TensorExpression<T,DIM,RANK,dataType> &tens) {
    (*this)(std::get<S>(indices)...) = tens(std::get<S>(indices)...);
  }

};

YAFEL_NAMESPACE_CLOSE

#endif
