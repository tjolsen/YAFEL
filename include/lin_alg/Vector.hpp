#ifndef _YAFEL_VECTOR_HPP
#define _YAFEL_VECTOR_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

class Vector {
  
private:
  static const unsigned default_capacity = 10;
  unsigned length;
  unsigned capacity;
  double *data;
  void resize();
  
public:
  Vector();
  Vector(unsigned len);
  Vector(unsigned len, double val);
  Vector(const Vector & src);
  ~Vector();
  Vector & operator=(const Vector & rhs);
  double & operator()(unsigned i) const;
  Vector& operator+=(const Vector & rhs);
  Vector operator+(const Vector & rhs) const;
  Vector operator-(const Vector & rhs) const;
  Vector & operator*=(double a);
  Vector operator*(double a) const; 
  void append(double val);
  double dot(const Vector & rhs) const;
  inline unsigned getLength() const {return length;}
  void print();
};

YAFEL_NAMESPACE_CLOSE

#endif
