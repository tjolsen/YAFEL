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
  value_type operator()(size_type i, size_type j) const {return _u(i,j)+_v(i,j)};
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
  value_type operator()(size_type i, size_type j) const {return _u(i,j)-_v(i,j)};
};



YAFEL_NAMESPACE_CLOSE

#endif
