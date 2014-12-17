#include <iostream>
#include "heap.hpp"

class Comparable {
  
private:
  int i;

public:
  Comparable(int I): i(I) {}
  bool operator<(const Comparable &rhs) {
    return (i < rhs.getVal());
  }
  int getVal() const { return i; }
  friend std::ostream& operator<<(std::ostream &os, const Comparable &c);
};

std::ostream& operator<<(std::ostream &os, const Comparable &c) {
  os << "Comparable(" << c.getVal() << ")";
  return os;
}

int main() {
  
  int N = 100;
  Heap<Comparable> h;
  
  for(int i=0; i<N; ++i) {
    std::cout << i << " ";
    Comparable c(N-i);
    h.insert(c);
  }
  std::cout << std::endl << std::endl;
  
  for(int i=0; i<N; ++i) {
    std::cout << h.extract() << std::endl;
  }
  
  return 0;
}
