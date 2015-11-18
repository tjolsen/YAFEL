#ifndef __YAFEL_TENSOREXPRESSION_HPP
#define __YAFEL_TENSOREXPRESSION_HPP

/*
 * TensorExpression:
 * Header file defining the expression templates for arbitrary rank/dimension tensors.
 * Classes use CRTP to allow the compiler to go to town optimizing nearly all function
 * calls away since all loop bounds and calls will be known at compile time.
 */


#include "yafel_globals.hpp"
#include "lin_alg/tensor/iterator_utils.hpp"
#include "lin_alg/tensor/reference_tensor_iterator.hpp"
#include "lin_alg/tensor/const_reference_tensor_iterator.hpp"
#include "lin_alg/tensor/reference_index_iterator.hpp"
#include "lin_alg/tensor/const_reference_index_iterator.hpp"

#include <tuple>
#include <type_traits>

YAFEL_NAMESPACE_OPEN

constexpr unsigned _tensorStorage(unsigned dim, unsigned rank) {
  return (rank==0) ? 1 : dim*_tensorStorage(dim,rank-1);
}

//--------------------------------------------------------------------------------
/*
 * TensorExpression base class. Repressents an arbitrary tensor of rank RANK and dimension DIM.
 * To be used via CRTP with classes representing tensor operations.
 */
//--------------------------------------------------------------------------------
template <typename T, unsigned DIM, unsigned RANK, typename dataType=double>
class TensorExpression {

public:
  using size_type = std::size_t;
  using value_type = dataType;

  size_type rank() const {return RANK;}
  size_type dim() const {return DIM;}

  template <typename ...Args>
  value_type operator()(Args ...args) const { return static_cast<const T&>(*this)(args...); }
  
  operator T&() {return static_cast<const T&>(*this);}
  operator T const&() const {return static_cast<const T&>(*this);}

  const_reference_tensor_iterator<TensorExpression<T,DIM,RANK,dataType>,DIM,RANK> begin() const {
    return const_reference_tensor_iterator<TensorExpression<T,DIM,RANK,dataType>,DIM,RANK>(*this);
  }

};


//--------------------------------------------------------------------------------
/*
 * Scalar scaling of a tensor expression. Scalar type T2 must be type-convertable to dataType
 */
//--------------------------------------------------------------------------------
template<typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType=double>
class TensorScaled : public TensorExpression<TensorScaled<T1,T2,DIM,RANK,dataType>, DIM, RANK, dataType> {
public:
  using value_type = typename TensorExpression<TensorScaled<T1,T2,DIM,RANK,dataType>,DIM,RANK,dataType>::value_type;
private:
  const T1 &_u;
  const dataType alpha;
  
public:
  TensorScaled(const TensorExpression<T1,DIM,RANK,dataType> &u, T2 a)
    : _u(u), alpha(dataType(a))
  {}
  
  template<typename ...Args>
  value_type operator()(Args ...args) const {return _u(args...)*alpha;}
};

template<typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType>
const TensorScaled<T1,T2,DIM,RANK,dataType>
operator*(const TensorExpression<T1,DIM,RANK,dataType> &u, T2 alpha) {
  static_assert(std::is_fundamental<T2>::value,
		"Error: TensorScaled currently only intended for primitive types");
  return TensorScaled<T1,T2,DIM,RANK,dataType>(u,alpha);
}

template<typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType>
const TensorScaled<T1,T2,DIM,RANK,dataType>
operator*(T2 alpha, const TensorExpression<T1,DIM,RANK,dataType> &u) {
  static_assert(std::is_fundamental<T2>::value,
		"Error: TensorScaled currently only intended for primitive types");
  return TensorScaled<T1,T2,DIM,RANK,dataType>(u,alpha);
}

//--------------------------------------------------------------------------------
/*
 * Class representing the sum of two TensorExpressions.
 * lhs (_u) and rhs (_v) must be of equal rank and dimension
 */
//--------------------------------------------------------------------------------
template <typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType=double>
class TensorSum : public TensorExpression<TensorSum<T1,T2,DIM,RANK,dataType>, DIM, RANK, dataType> {
public:
  using value_type = typename TensorExpression<TensorSum<T1,T2,DIM,RANK,dataType>,DIM,RANK,dataType>::value_type;
private:
  const T1 &_u;
  const T2 &_v;

public:

  TensorSum(const TensorExpression<T1,DIM,RANK,dataType> &u,
	    const TensorExpression<T2,DIM,RANK,dataType> &v)
    : _u(u), _v(v)
  {}
  
  template<typename ...Args>
  value_type operator()(Args ...args) const {return _u(args...)+_v(args...);}

};

template<typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType>
const TensorSum<T1,T2,DIM,RANK,dataType>
operator+(const TensorExpression<T1,DIM,RANK,dataType> &u,
	  const TensorExpression<T2,DIM,RANK,dataType> &v) {
  return TensorSum<T1,T2,DIM,RANK,dataType>(u,v);
}

//--------------------------------------------------------------------------------
/*
 * Class representing the difference between two TensorExpressions.
 * lhs (_u) and rhs (_v) must be of equal rank and dimension
 */
//--------------------------------------------------------------------------------
template <typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType=double>
class TensorDifference : 
  public TensorExpression<TensorDifference<T1,T2,DIM,RANK,dataType>, DIM, RANK, dataType> 
{
public:
  using value_type = typename TensorExpression<TensorDifference<T1,T2,DIM,RANK,dataType>,DIM,RANK,dataType>::value_type;
private:
  const T1 &_u;
  const T2 &_v;

public:

  TensorDifference(const TensorExpression<T1,DIM,RANK,dataType> &u,
		   const TensorExpression<T2,DIM,RANK,dataType> &v)
    : _u(u), _v(v)
  {}
  
  template<typename ...Args>
  value_type operator()(Args ...args) const {return _u(args...)-_v(args...);}

};

template<typename T1, typename T2, unsigned DIM, unsigned RANK, typename dataType=double>
const TensorDifference<T1,T2,DIM,RANK,dataType>
operator-(const TensorExpression<T1,DIM,RANK,dataType> &u,
	  const TensorExpression<T2,DIM,RANK,dataType> &v) {
  return TensorDifference<T1,T2,DIM,RANK,dataType>(u,v);
}


//--------------------------------------------------------------------------------
/*
 * Tensor outer product. Result is a new tensor expression of dimension DIM, rank R1+R2.
 */
//--------------------------------------------------------------------------------
template <typename T1, typename T2, unsigned DIM, unsigned R1, unsigned R2, typename dataType=double>
class OuterProduct : public TensorExpression<OuterProduct<T1,T2,DIM,R1,R2,dataType>,DIM,R1+R2,dataType>
{
public:
  using value_type = typename TensorExpression<OuterProduct<T1,T2,DIM,R1,R2,dataType>,DIM,R1+R2,dataType>::value_type;
private:
  const T1 & _u;
  const T2 & _v;
  
  template <int ...S, typename ...Args>
  value_type lhs(seq<S...>, const std::tuple<Args...> &params) const {
    return _u(std::get<S>(params)...);
  }
  
  template <int ...S, typename ...Args>
  value_type rhs(seq<S...>, const std::tuple<Args...> &params) const {
    return _v(std::get<S+R1>(params)...);
  }
  
public:
  OuterProduct(const TensorExpression<T1,DIM,R1,dataType> &u, 
	       const TensorExpression<T2,DIM,R2,dataType> &v) 
  : _u(u), _v(v)
  {}
  
  template<typename ...Args>
  value_type operator()(Args ...args) const {
    static_assert(sizeof...(args) == R1+R2, 
		  "OuterProduct::operator() error wrong number of arguments");

    return lhs(typename gens<R1>::type(), std::forward_as_tuple(args...))
      *rhs(typename gens<R2>::type(), std::forward_as_tuple(args...));
  }
  
};

template<typename T1, typename T2, unsigned DIM, unsigned R1, unsigned R2, typename dataType=double>
const OuterProduct<T1,T2,DIM,R1,R2,dataType>
otimes(const TensorExpression<T1, DIM, R1, dataType> &u,
       const TensorExpression<T2, DIM, R2, dataType> &v) 
{
  return OuterProduct<T1, T2, DIM, R1, R2, dataType>(u,v);
}


//--------------------------------------------------------------------------------
/*
 * Class representing the inner-product contraction of 2 tensors.
 * The NCONTRACT template parameter specifies how many indices are to be contracted out
 *
 * Need to figure out way to specialize to return a dataType rather than a TensorExpression
 * for operations resulting in a rank 0 tensor.
 *
 * Example operations would be: 
 *                             c_i = A_{ij}b_{j}      <--- NCONTRACT = 1
 *                             A_{ij} = C_{ijkl}B_{kl}    <--- NCONTRACT = 2
 */
//--------------------------------------------------------------------------------
template <typename T1, typename T2, unsigned DIM, unsigned R1, unsigned R2, 
	  unsigned NCONTRACT, typename dataType=double>
class TensorContraction :
  public TensorExpression<TensorContraction<T1,T2,DIM,R1,R2,NCONTRACT,dataType>,DIM,R1+R2-2*NCONTRACT, dataType>
{
public:
  using value_type = typename TensorExpression<TensorContraction<T1,T2,DIM,R1,R2,NCONTRACT,dataType>,
					       DIM, R1+R2-2*NCONTRACT, dataType>::value_type;

private:
  const T1 &_u;
  const T2 &_v;

public:
  
  TensorContraction(const TensorExpression<T1,DIM,R1,dataType> &u,
		    const TensorExpression<T2,DIM,R2,dataType> &v)
    : _u(u), _v(v)
  {}


};



YAFEL_NAMESPACE_CLOSE

#endif
