#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/generic_tensor_iterator.hpp"
#include "lin_alg/tensor/generic_index_iterator.hpp"

#include <iostream>

using namespace yafel;

// Test that generic_tensor_iterator is hitting the right number of locations
bool test_1() {

  generic_tensor_iterator<3,4> GI;
  std::size_t count=0;
  
  for(; !GI.end(); GI.next()) {
    ++count;
  }
  
  return count == 3*3*3*3;
}


// Test that generic_index_iterator hits correct number of locations
bool test_2() {
  generic_index_iterator<3,4,0,1> GI(0,0,0,0);
  std::size_t count=0;
  
  for(; !GI.end(); GI.next()) {
    ++count;
  }
  return count==3*3;
}

// Test construction of Tensor object, assignment of values, then sum of values
bool test_3() {

  Tensor<3,2> A;
  
  int a = 1;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = a;
    ++a;
  }

  int sum=0;
  for(auto it = A.begin(); !it.end(); it.next()) {
    sum += *it;
  }
  return sum == 45;
}

// Test tensor addition
bool test_4() {
  Tensor<3,2,int> A,B;

  int a=1, b=2;

  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = a;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = b;
  }

  Tensor<3,2,int> C = A+B;
  bool good = true;
  for(auto it=C.begin(); !it.end(); it.next()) {
    good = good && (*it==(a+b));
  }
  
  return good;
}

// Test tensor subtraction
bool test_5() {
  Tensor<3,5,int> A,B;

  int a=1, b=2;

  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = a;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = b;
  }

  Tensor<3,5,int> C = A-B;
  bool good = true;
  for(auto it=C.begin(); !it.end(); it.next()) {
    good = good && (*it==(a-b));
  }
  
  return good;
}

// Test tensor scaling
bool test_6() {

  Tensor<4,4,int> A,B;

  // put together some contrived calculation that couldn't be accidental
  int a = 4, b = -3;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = 2;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = 1;
  }
  
  Tensor<4,4,int> C = a*A + b*B;

  bool good = true;
  for(auto it=C.begin(); !it.end(); it.next()) {
    good = good && (*it == (2*a + b));
  }

  return good;
}

// test outer product
bool test_7() {
  
  Tensor<3,3,double> A;
  Tensor<3,2,double> B;
  
  double a=2, b=3;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it=a;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it=b;
  }

  bool good = true;
  Tensor<3,5> C = otimes(A,B);
  for(auto it=C.begin(); !it.end(); it.next()) {
    good = good && (*it == a*b);
  }

  return good;
}

// Test tensor-contraction of rank 2 tensors (ie: matmul)
bool test_8() {
  
  Tensor<3,2> A,B;
  
  for(std::size_t i=0; i<3; ++i) {
    for(std::size_t j=0; j<3; ++j)  {
      A(i,j) = 1.0;
      B(i,j) = 2.0;
    }
  }
  
  Tensor<3,2> C = contract<1>(A,B);
  
  
  return C(0,0) == 6.0;
}


int main() {

  int retval = 0;
  
  if(!test_1()) {
    retval |= 1<<0;
    std::cout << "Failed test_1" << "\n";
  }
  if(!test_2()) {
    retval |= 1<<1;
    std::cout << "Failed test_2" << "\n";
  }
  if(!test_3()) {
    retval |= 1<<2;
    std::cout << "Failed test_3" << "\n";
  }
  if(!test_4()) {
    retval |= 1<<3;
    std::cout << "Failed test_4" << "\n";
  }
  if(!test_5()) {
    retval |= 1<<4;
    std::cout << "Failed test_5" << "\n";
  }
  if(!test_6()) {
    retval |= 1<<5;
    std::cout << "Failed test_6" << "\n";
  }
  if(!test_7()) {
    retval |= 1<<6;
    std::cout << "Failed test_7" << "\n";
  }
  if(!test_8()) {
    retval |= 1<<7;
    std::cout << "Failed test_8" << "\n";
  }

  return retval;
}
