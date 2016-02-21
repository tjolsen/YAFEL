#ifndef __YAFEL_TENSOR_HPP
#define __YAFEL_TENSOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"
#include <initializer_list>

YAFEL_NAMESPACE_OPEN

template<unsigned DIM, unsigned RANK, typename dataType=double>
class Tensor : public TensorExpression< Tensor<DIM,RANK,dataType>, DIM, RANK, dataType> {
public:
  using value_type = typename TensorExpression<Tensor<DIM,RANK,dataType>,DIM,RANK,dataType>::value_type;
  using size_type = typename TensorExpression<Tensor<DIM,RANK,dataType>,DIM,RANK,dataType>::size_type;
  
  template<typename ...Args>
  value_type & operator()(Args ...args) {return _data[index(args...)];}

  template<typename ...Args>
  value_type operator()(Args ...args) const {return _data[index(args...)];}
  
  //get a reference_tensor_iterator to the beginning
  reference_tensor_iterator<Tensor<DIM,RANK,dataType>,DIM,RANK> begin() {
    return reference_tensor_iterator<Tensor<DIM,RANK,dataType>,DIM,RANK>(*this);
  }


  //initializer list construction for rank 1 tensors
  Tensor(std::initializer_list<dataType> il) {
    static_assert(RANK==1, "Tensor: Tried to use initializer list ctor for RANK \\neq 1");

#ifndef _OPTIMIZED
    if(il.size() > DIM) {
      throw(std::length_error("Tensor: Error: Initializer list ctor size() \\neq DIM"));
    }
#endif

    size_type i=0;
    for(auto it=il.begin(); it<il.end(); ++i, ++it) {
      (*this)(i) = *it;
    }

  }

  
  //initialize a zero-tensor
  Tensor() {

    for(auto ti = this->begin(); !ti.end(); ti.next()) {
      *ti = value_type(0);
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


  /*
   * Update operators +=, -=
   */
  template<typename T>
  Tensor<DIM,RANK,dataType> & operator+=(const TensorExpression<T,DIM,RANK,dataType> &rhs) {
    
    auto lhsit = this->begin();
    auto rhsit = rhs.begin();
    for(; !lhsit.end(); lhsit.next(), rhsit.next()) {
      *lhsit += *rhsit;
    }
    
    return *this;
  }

  template<typename T>
  Tensor<DIM,RANK,dataType> & operator-=(const TensorExpression<T,DIM,RANK,dataType> &rhs) {
    
    auto lhsit = this->begin();
    auto rhsit = rhs.begin();
    for(; !lhsit.end(); lhsit.next(), rhsit.next()) {
      *lhsit -= *rhsit;
    }

    return *this;
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

};


YAFEL_NAMESPACE_CLOSE

#endif
