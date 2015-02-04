#ifndef _HEAP_HPP
#define _HEAP_HPP

#include <vector>
#include <iostream>

//Template class for a binary heap


template<typename T>
class Heap {
  
private:
  std::vector<T> array;

  void BubbleUp(unsigned node);
  void BubbleDown(unsigned node);
  inline unsigned childL(unsigned node) {return node*2 + 1;}
  inline unsigned childR(unsigned node) {return node*2 + 2;}
  inline unsigned parent(unsigned node) {return (node-1)/2;}
  
public:
  Heap();
  void insert(T item);
  T peek();
  T extract();
  void print();
  inline unsigned size() const {return array.size();}
};

#include "heap_implementation.hpp"

#endif
