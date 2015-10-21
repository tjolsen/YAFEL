#ifndef __YAFEL_VECTOR_HPP
#define __YAFEL_VECTOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/VectorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class Vector: public VectorExpression<Vector<dataType>, dataType> {
  
private:
  typedef typename VectorExpression<Vector<dataType>,dataType>::container_type container_type;
  typedef typename VectorExpression<Vector<dataType>,dataType>::value_type value_type;
  typedef typename VectorExpression<Vector<dataType>,dataType>::size_type size_type;
  typedef typename VectorExpression<Vector<dataType>,dataType>::reference reference;
  container_type _data;


public:
  reference operator()(size_type i) {return _data[i];}
  value_type operator()(size_type i) const {return _data[i];}
  size_type size() const {return _data.size();}
  
  /*
   * Constructors
   */
  // Construct Vector of length N
  Vector(size_type N) : _data(N) {}

  // Construct Vector of length N and fill with value
  Vector(size_type N, value_type val) : _data(N, val) {}
  
  // Construct from VectorExpression<T,dataType>
  template<typename T>
  Vector(const VectorExpression<T,dataType> & v) {
    _data.resize(v.size());
    for(size_type i=0; i<v.size(); ++i) {
      _data[i] = v(i);
    }
  }

};


YAFEL_NAMESPACE_CLOSE

#endif
