#ifndef _SPATIALFUNCTION_HPP
#define _SPATIALFUNCTION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T>
class SpatialFunction {

private:
  T (*func)(const Vector &x);

public:
  SpatialFunction(T (*fp)(const Vector &x)) : func(fp) {}
  T operator()(const Vector &x) const {return func(x);}
};

YAFEL_NAMESPACE_CLOSE

#endif
