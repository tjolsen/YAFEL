#include <iostream>
#include "heap.hpp"


int main() {
  
  int N = 10;
  Heap<int> h;
  
  for(int i=0; i<N; ++i) {
    std::cout << i << " ";
    h.insert(10-i);
  }
  std::cout << std::endl;
  
  return 0;
}
