#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/generic_tensor_iterator.hpp"
#include "lin_alg/tensor/generic_index_iterator.hpp"

#include <iostream>

using namespace yafel;

bool test_1() {

  generic_tensor_iterator<3,2> GI;
  std::size_t count=0;
  
  for(; !GI.end(); GI.next()) {
    ++count;
  }
  
  return count == 3*3;
}



int main() {

  int retval = 0;
  
  if(!test_1()) {
    retval |= 1<<0;
    std::cout << "Failed test_1" << "\n";
  }

  return retval;
}
