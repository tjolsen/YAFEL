#ifndef __YAFEL_MATRIXEXPRESSION_HPP
#define __YAFEL_MATRIXEXPRESSION_HPP

#include "yafel_globals.hpp"

#include <vector>
#include <cassert>

YAFEL_NAMESPACE_OPEN

/*
 * MatrixExpression: Base class utilizing CRTP to define efficient matrix operations at compile-time.
 * This class cannot be instantiated directly.
 */
template<typename T, typename dataType=double>
class MatrixExpression {
public:
  typedef std::vector<dataType> container_type;
  typedef typename container_type::value_type value_type;
  typedef typename container_type::size_type size_type;
  typedef typename container_type::reference reference;

  size_type rows() const {return static_cast<T const&>(*this).rows();}
  size_type cols() const {return static_cast<T const&>(*this).cols();}
  value_type operator()(size_type i, size_type j) const {return static_cast<const T&>(*this)(i,j);}
  
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() const {return static_cast<const T&>(*this);}
};



/*
 * Class definitions representing specific operations(+, -, scalar*, ^T, ...)
 */
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType=double>
class MatrixSum : public MatrixExpression<MatrixSum<T1,T2,dataType>, dataType> {

private:
  const T1 &_u;
  const T2 &_v;

public:
  typedef typename MatrixExpression<MatrixSum<T1,T2,dataType>, dataType>::value_type value_type;
  typedef typename MatrixExpression<MatrixSum<T1,T2,dataType>, dataType>::size_type size_type;

  /*
   * Constructor
   */
  MatrixSum(const MatrixExpression<T1,dataType> &u,
	    const MatrixExpression<T2,dataType> &v) : _u(u), _v(v) {
    assert(u.rows() == v.rows() && "MatrixSum: rows mismatch");
    assert(u.cols() == v.cols() && "MatrixSum: cols mismatch");
  }
  
  size_type rows() const {return _u.rows();}
  size_type cols() const {return _u.cols();}
  value_type operator()(size_type i, size_type j) const {return _u(i,j)+_v(i,j);}
};

//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType=double>
class MatrixDifference : public MatrixExpression<MatrixDifference<T1,T2,dataType>, dataType> {

private:
  const T1 &_u;
  const T2 &_v;

public:
  typedef typename MatrixExpression<MatrixDifference<T1,T2,dataType>, dataType>::value_type value_type;
  typedef typename MatrixExpression<MatrixDifference<T1,T2,dataType>, dataType>::size_type size_type;

  /*
   * Constructor
   */
  MatrixDifference(const MatrixExpression<T1,dataType> &u,
		   const MatrixExpression<T2,dataType> &v) : _u(u), _v(v) {
    assert(u.rows() == v.rows() && "MatrixDifference: rows mismatch");
    assert(u.cols() == v.cols() && "MatrixDifference: cols mismatch");
  }
  
  size_type rows() const {return _u.rows();}
  size_type cols() const {return _u.cols();}
  value_type operator()(size_type i, size_type j) const {return _u(i,j)-_v(i,j);}
};

//----------------------------------------------------------------------------------------
template<typename T1, typename dataType=double>
class MatrixScaled : public MatrixExpression<MatrixScaled<T1,dataType>, dataType> {

private:
  const dataType scalar;
  const T1 &_u;

public:
  typedef typename MatrixExpression<MatrixScaled<T1,dataType>, dataType>::value_type value_type;
  typedef typename MatrixExpression<MatrixScaled<T1,dataType>, dataType>::size_type size_type;

  /*
   * Constructor
   */
  template<typename T2>
  MatrixScaled(const MatrixExpression<T1,dataType> &u, T2 alpha) : scalar((dataType)alpha), _u(u) {}

  size_type rows() const {return _u.rows();}
  size_type cols() const {return _u.cols();}
  value_type operator()(size_type i, size_type j) const {return scalar*_u(i,j);}
};

//----------------------------------------------------------------------------------------
template<typename T1, typename dataType=double>
class MatrixTranspose : public MatrixExpression<MatrixTranspose<T1,dataType>, dataType> {
private:
  const T1 &_u;
  
public:
  typedef typename MatrixExpression<MatrixTranspose<T1,dataType>,dataType>::value_type value_type;
  typedef typename MatrixExpression<MatrixTranspose<T1,dataType>,dataType>::size_type size_type;

  /*
   * Constructor
   */
  MatrixTranspose(const MatrixExpression<T1,dataType> &u) : _u(u) {}

  size_type rows() const {return _u.rows();}
  size_type cols() const {return _u.cols();}
  value_type operator()(size_type i, size_type j) const {return u(j,i);}
};


//----------------------------------------------------------------------------------------
/*
 * Operator overloading: Define operators (+, -, *, ...) for MatrixExpressions
 */
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
MatrixSum<T1,T2,dataType> operator+(const MatrixExpression<T1,dataType> &lhs,
				    const MatrixExpression<T2,dataType> &rhs) {
  return MatrixSum<T1,T2,dataType>(lhs,rhs);
}
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
MatrixDifference<T1,T2,dataType> operator-(const MatrixExpression<T1,dataType> &lhs,
					   const MatrixExpression<T2,dataType> &rhs) {
  return MatrixDifference<T1,T2,dataType>(lhs,rhs);
}
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
MatrixScaled<T1,dataType> operator*(const MatrixExpression<T1,dataType> &v, T2 a) {
  return MatrixScaled<T1,dataType>(v,a);
}
template<typename T1, typename T2, typename dataType>
MatrixScaled<T1,dataType> operator*(T2 a, const MatrixExpression<T1,dataType> &v) {
  return MatrixScaled<T1,dataType>(v,a);
}

//----------------------------------------------------------------------------------------
/*
 * Some relational operators that may be userful
 */
//----------------------------------------------------------------------------------------
template<typename T1, typename T2, typename dataType>
bool operator==(const MatrixExpression<T1,dataType> &lhs,
		const MatrixExpression<T2,dataType> &rhs) {
  
  if(lhs.rows()!=rhs.rows() || lhs.cols()!=rhs.cols()) {
    return false;
  }

  for(typename MatrixExpression<T1,dataType>::size_type i=0; i<lhs.rows(); ++i) {
    for(typename MatrixExpression<T1,dataType>::size_type j=0; j<lhs.cols(); ++j) {
      if(lhs(i,j) != rhs(i,j)) {
	return false;
      }
    }
  }
  
  return true;
}

YAFEL_NAMESPACE_CLOSE

#endif
