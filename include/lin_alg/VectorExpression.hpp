#ifndef __YAFEL_VECTOREXPRESSION_HPP
#define __YAFEL_VECTOREXPRESSION_HPP

#include "yafel_globals.hpp"

#include <vector>
#include <cassert>

YAFEL_NAMESPACE_OPEN



/*
 * VectorExpression: Base class utilizing CRTP to define efficient vector operations at compile-time.
 * This class cannot be instantiated directly.
 */
template<typename T, typename dataType=double>
class VectorExpression {
public:
  typedef std::vector<dataType> container_type;
  typedef typename container_type::value_type value_type;
  typedef typename container_type::size_type size_type;
  typedef typename container_type::reference reference;

  size_type size() const { return static_cast<T const&>(*this).size(); }
  value_type operator[](size_type i) const { return static_cast<const T&>(*this)[i]; }
  
  operator T&() { return static_cast<T&>(*this); }
  operator T const&() const { return static_cast<const T&>(*this); }
};


/*
 * Class definition representing specific operations (+, -, *,...)
 */

//-------------------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType=double>
class VectorDifference : public VectorExpression<VectorDifference<T1,T2,dataType>, dataType> {

private:
  const T1 &_u;
  const T2 &_v;

public:
  typedef typename VectorExpression<VectorDifference<T1,T2,dataType>,dataType>::value_type value_type;
  typedef typename VectorExpression<VectorDifference<T1,T2,dataType>,dataType>::size_type size_type;

  /*
   * Constructor
   */
  VectorDifference(const VectorExpression<T1> &u, const VectorExpression<T2> &v) : _u(u), _v(v) {
    assert(u.size() == v.size());
  }
  
  value_type operator[](size_type i) const {return _u[i] - _v[i];}
  size_type size() const {return _u.size();}
};


//-------------------------------------------------------------------------------------------------
/*
 * VectorSum: Class to represent the element-wise addition of two VectorExpressions
 */
template<typename T1, typename T2, typename dataType=double>
class VectorSum : public VectorExpression<VectorSum<T1,T2,dataType>, dataType> {

private:
  const T1 &_u;
  const T2 &_v;

public:
  typedef typename VectorExpression<VectorSum<T1,T2,dataType>,dataType>::value_type value_type;
  typedef typename VectorExpression<VectorSum<T1,T2,dataType>,dataType>::size_type size_type;

  /*
   * Constructor
   */
  VectorSum(const VectorExpression<T1,dataType> &u, const VectorExpression<T2,dataType> &v) : _u(u), _v(v) {
    assert(u.size() == v.size());
  }
  
  /*
   * Member functions
   */
  value_type operator[](size_type i) const {return _u[i] + _v[i];}
  size_type size() const {return _v.size();}
};


//-------------------------------------------------------------------------------------------------
/*
 * VectorScaled: Class to represent the multiplicative scaling of a VectorExpression by a scalar
 */
template<typename T1, typename dataType=double>
class VectorScaled : public VectorExpression<VectorScaled<T1,dataType>, dataType> {
private:
  const dataType _scalar;
  const T1 &_v;

public:
  typedef typename VectorExpression<VectorScaled<T1,dataType>,dataType>::value_type value_type;
  typedef typename VectorExpression<VectorScaled<T1,dataType>,dataType>::size_type size_type;
  
  /*
   * Constructor
   */
  template<typename T2>
  VectorScaled(const VectorExpression<T1,dataType> &v, T2 alpha) : _scalar((dataType)alpha), _v(v) {}
  
  /*
   * Member functions
   */
  value_type operator[](size_type i) const {return _scalar*_v[i];}
  size_type size() const {return _v.size();}
};


//-------------------------------------------------------------------------
/*
 * Operator overloading: Define operators (+,-,*,...) for VectorExpressions
 */
//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
VectorSum<T1,T2,dataType> operator+(const VectorExpression<T1,dataType> &lhs,
				    const VectorExpression<T2,dataType> &rhs) {
  return VectorSum<T1,T2,dataType>(lhs,rhs);
}

//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
VectorDifference<T1,T2,dataType> operator-(const VectorExpression<T1,dataType> &lhs,
				    const VectorExpression<T2,dataType> &rhs) {
  return VectorDifference<T1,T2,dataType>(lhs,rhs);
}

//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
VectorScaled<T1,dataType> operator*(const VectorExpression<T1,dataType> &v, T2 a) {
  return VectorScaled<T1,dataType>(v, a);
}
template<typename T1, typename T2, typename dataType>
VectorScaled<T1,dataType> operator*(T2 a, const VectorExpression<T1,dataType> &v) {
  return VectorScaled<T1,dataType>(v, a);
}


//-------------------------------------------------------------------------
/*
 * Some relational operators (may or may not be useful. At least useful for testing.
 */
//-------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
bool operator==(const VectorExpression<T1,dataType> &lhs,
		const VectorExpression<T2,dataType> &rhs) {
  if(lhs.size() != rhs.size()) {
    return false;
  }
  
  for(typename VectorExpression<T1,dataType>::size_type i=0; i<lhs.size(); ++i) {
    if(lhs[i] != rhs[i]) {
      return false;
    }
  }
  
  return true;
}

YAFEL_NAMESPACE_CLOSE

#endif
