#ifndef _YAFEL_VECTOR_HPP
#define _YAFEL_VECTOR_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

class Vector {
  
private:
  static const int default_capacity = 10;
  int length;
  int capacity;
  double *data;
  void resize();
  
public:
  Vector();
  Vector(int len);
  Vector(int len, double val);
  Vector(const Vector & src);
  ~Vector();
  Vector & operator=(const Vector & rhs);
  double & operator()(int i) const;
  Vector& operator+=(const Vector & rhs);
  Vector operator+(const Vector & rhs) const;
  Vector operator-(const Vector & rhs) const;
  Vector & operator*=(double a);
  Vector operator*(double a) const; 
  void append(double val);
  double dot(const Vector & rhs) const;
  inline int getLength() const {return length;}
  void print();
};

YAFEL_NAMESPACE_CLOSE

#endif
