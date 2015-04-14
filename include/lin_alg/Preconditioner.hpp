#ifndef _YAFEL_PRECONDITIONER_HPP
#define _YAFEL_PRECONDITIONER_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

class Preconditioner {


public:
  //all preconditioners must implement a (hopefully fast) linear solve.
  virtual Vector MinvV(const Vector &rhs) const = 0;
  

};

YAFEL_NAMESPACE_CLOSE

#endif
