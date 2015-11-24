#ifndef __YAFEL_DUALNUMBER_HPP
#define __YAFEL_DUALNUMBER_HPP

#include "yafel_globals.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

template<typename T>
class DualNumber {
public:
  T first;
  T second;
  DualNumber():DualNumber(0,0) {}
  DualNumber(T v1, T v2): first(v1), second(v2) {}

  // arithmetic operator overloading (+, -, *, /)
  DualNumber<T> operator+(const DualNumber<T> &rhs) const {
    return DualNumber<T>(first+rhs.first, second+rhs.second);
  }
  DualNumber<T> operator+(double rhs) const {
    return DualNumber<T>(first+rhs, second);
  }

  DualNumber<T> operator-(const DualNumber<T> & rhs) const {
    return DualNumber<T>(first-rhs.first, second-rhs.second);
  }
  DualNumber<T> operator-(double rhs) const {
    return DualNumber<T>(first-rhs, second);
  }
  
  DualNumber<T> operator*(const DualNumber<T> & rhs) const {
    return DualNumber<T>(first*rhs.first, second*rhs.first + first*rhs.second);
  }
  DualNumber<T> operator*(double rhs) const {
    return DualNumber<T>(first*rhs, second*rhs);
  }

  DualNumber<T> operator/(const DualNumber<T> &rhs) const {
    return DualNumber<T>(first/rhs.first, 
			 (first*rhs.second - second*rhs.first)/(rhs.first*rhs.first));
  }
  DualNumber<T> operator/(double rhs) const {
    return DualNumber<T>(first/rhs, second/rhs);
  }

  
  // unary operator-()
  DualNumber<T> operator-() const {
    return DualNumber<T>(-first, -second);
  }
};

// Overload operators for primitive types with lhs/rhs order reversed
template<typename T>
DualNumber<T> operator+(double lhs, DualNumber<T> rhs) {
  return DualNumber<T>(lhs+rhs.first, rhs.second);
}
template<typename T>
DualNumber<T> operator-(double lhs, DualNumber<T> rhs) {
  return DualNumber<T>(lhs-rhs.first, rhs.second);
}
template<typename T>
DualNumber<T> operator*(double lhs, DualNumber<T> rhs) {
  return DualNumber<T>(lhs*rhs.first, lhs*rhs.second);
}
template<typename T>
DualNumber<T> operator/(double lhs, DualNumber<T> rhs) {
  return DualNumber<T>(lhs/rhs.first, -lhs*rhs.second/(rhs.first*rhs.first));
}

// More Useful Functions
template<typename T>
DualNumber<T> sin(DualNumber<T> x) {
  return DualNumber<T>(sin(x.first), x.second*cos(x.first));
}

template<typename T>
DualNumber<T> cos(DualNumber<T> x) {
  return DualNumber<T>(cos(x.first), -x.second*sin(x.first));
}

template<typename T>
DualNumber<T> exp(DualNumber<T> x) {
  return DualNumber<T>(exp(x.first), x.second*exp(x.first));
}


YAFEL_NAMESPACE_CLOSE

#endif
