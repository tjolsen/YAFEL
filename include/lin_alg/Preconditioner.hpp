#ifndef _YAFEL_PRECONDITIONER_HPP
#define _YAFEL_PRECONDITIONER_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T, typename dataType>
class Preconditioner {
public:
  operator T&() {return static_cast<T&>(*this);}
  operator T const&() {return static_cast<T const&>(*this);}
  
  //all preconditioners must implement a (hopefully fast) linear solve.
  Vector<dataType> MinvV(const Vector<double> &rhs) const {
    return static_cast<const T&>(*this).MinvV(rhs);
  }

};

YAFEL_NAMESPACE_CLOSE

#endif
