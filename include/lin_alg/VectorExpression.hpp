#ifndef _YAFEL_VECTOREXPRESSION_HPP
#define _YAFEL_VECTOREXPRESSION_HPP

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
  value_type operator()(size_type i) const { return static_cast<const T&>(*this)(i); }
  
  operator T&() { return static_cast<T&>(*this); }
  operator T const&() const { return static_cast<const T&>(*this); }

  //dot product
  template<typename T2>
  value_type dot(const VectorExpression<T2,dataType> &rhs) const {
#ifndef _OPTIMIZED
    assert(this->size() == rhs.size() && "VectorExpression::dot dimension mismatch");
#endif
    
    value_type retval = value_type(0);
    for(size_type i=0; i<size(); ++i) {
      retval += (*this)(i)*rhs(i);
    }
    
    return retval;
  }
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
  VectorDifference(const VectorExpression<T1,dataType> &u, 
		   const VectorExpression<T2,dataType> &v) : _u(u), _v(v) {
    assert(u.size() == v.size());
  }
  
  value_type operator()(size_type i) const {return _u(i) - _v(i);}
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
  value_type operator()(size_type i) const {return _u(i) + _v(i);}
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
  value_type operator()(size_type i) const {return _scalar*_v(i);}
  size_type size() const {return _v.size();}
};



YAFEL_NAMESPACE_CLOSE

#endif
