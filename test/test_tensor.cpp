#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/generic_tensor_iterator.hpp"
#include "lin_alg/tensor/generic_index_iterator.hpp"
#include "lin_alg/tensor/tensor_specializations.hpp"

#include <iostream>
#include <typeinfo>

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

  auto C = A+B; //Tensor<3,2,int>
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

  auto C = A-B; //Tensor<3,5,int>
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
  
  Tensor<3,3,int> A;
  Tensor<3,2,int> B;
  
  double a=2, b=3;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it=a;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it=b;
  }

  bool good = true;
  Tensor<3,5,int> C = otimes(A,B);
  for(auto it=C.begin(); !it.end(); it.next()) {
    good = good && (*it == a*b);
  }

  return good;
}

// Test tensor-contraction of rank 2 tensors (ie: matmul)
bool test_8() {
  
  Tensor<3,2,int> A,B;
  
  int counter = 3;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = counter;
    counter = (7*counter)%13;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = counter;
    counter = (7*counter)%13;
  }
  
  Tensor<3,2,int> C = contract<1>(A,B);
  
  Tensor<3,2,int> D;
  bool good = true;
  for(std::size_t i=0; i<3; ++i) {
    for(std::size_t j=0; j<3; ++j) {
      D(i,j) = 0;
      for(std::size_t k=0; k<3; ++k) {
	D(i,j) += A(i,k)*B(k,j);
      }
      good = good && C(i,j)==D(i,j);
    }
  }

  return good;
}


// Test tensor-contraction of rank 4 tensors
bool test_9() {
  
  Tensor<3,4,int> A,B;
  int a = 7;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }
  

  Tensor<3,4,int> C = contract<2>(A,B);
  Tensor<3,4,int> D;
  bool good = true;
  for(std::size_t i=0; i<3; ++i) {
    for(std::size_t j=0; j<3; ++j) {
      for(std::size_t k=0; k<3; ++k) {
	for(std::size_t l=0; l<3; ++l) {
	  D(i,j,k,l) = 0;
	  for(std::size_t m=0; m<3; ++m) {
	    for(std::size_t n=0; n<3; ++n) {
	      D(i,j,k,l) += A(i,j,m,n)*B(m,n,k,l);
	    }
	  }
	  good = good && C(i,j,k,l) == D(i,j,k,l);
	}
      }
    }
  }

  return good;
}

// Test tensor contraction of rank 2 and rank 1 tensors
bool test_10() {
  Tensor<3,2,int> A;
  Tensor<3,1,int> B;

  int a = 7;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }

  Tensor<3,1,int> C = contract<1>(A,B);
  Tensor<3,1,int> D;
  bool good = true;
  for(std::size_t i=0; i<3; ++i) {
    D(i) = 0;
    for(std::size_t j=0; j<3; ++j) {
      D(i) += A(i,j)*B(j);
    }
    good = good && C(i)==D(i);
  }
  
  return good;
}


// Test tensor contraction of rank 2 and rank 1 tensor expressions
bool test_11() {
  Tensor<3,2,int> A, B;
  Tensor<3,1,int> u,v;

  int a = 7;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }
  for(auto it=u.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }
  for(auto it=v.begin(); !it.end(); it.next()) {
    *it = a;
    a = (17*a)%31;
  }

  Tensor<3,1,int> C = contract<1>(3*A+B, u-v*2);

  Tensor<3,1,int> D;
  bool good = true;
  for(std::size_t i=0; i<3; ++i) {
    D(i) = 0;
    for(std::size_t j=0; j<3; ++j) {
      D(i) += (3*A(i,j) + B(i,j))*(u(j) - v(j)*2);
    }
    good = good && C(i)==D(i);
  }
  
  return good;
}


// Test tensor-contraction of rank 2 tensors (ie: matmul) using operator*
bool test_12() {
  
  Tensor<3,2,int> A,B;
  
  int counter = 3;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = counter;
    counter = (7*counter)%13;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = counter;
    counter = (7*counter)%13;
  }
  
  Tensor<3,2,int> C = A*B;
  
  Tensor<3,2,int> D;
  bool good = true;
  for(std::size_t i=0; i<3; ++i) {
    for(std::size_t j=0; j<3; ++j) {
      D(i,j) = 0;
      for(std::size_t k=0; k<3; ++k) {
	D(i,j) += A(i,k)*B(k,j);
      }
      good = good && C(i,j)==D(i,j);
    }
  }

  return good;
}

// Test associativity of tensor multiplication: (a*A)*B == a*(A*B)
bool test_13() {
  
  Tensor<3,2,int> A,B;

  int a = -5;
  
  int counter = 3;
  for(auto it=A.begin(); !it.end(); it.next()) {
    *it = counter;
    counter = (7*counter)%13;
  }
  for(auto it=B.begin(); !it.end(); it.next()) {
    *it = counter;
    counter = (7*counter)%13;
  }
  
  Tensor<3,2,int> C = (a*A)*B;
  
  Tensor<3,2,int> D = a*(A*B);

  bool good = true;
  for(std::size_t i=0; i<3; ++i) {
    for(std::size_t j=0; j<3; ++j) {
      good = good && C(i,j)==D(i,j);
    }
  }

  return good;
}


// test full contraction of rank 1 tensors (dot product)
bool test_14() {
  Tensor<3,1,int> A, B;
  auto ait=A.begin();
  auto bit=B.begin();
  for(; !ait.end(); ait.next(), bit.next()) {
    *ait = 1;
    *bit = 1;
  }

  return contract<1>(A,B)==3;
}

// test full contraction of rank 4 tensors (dot product)
bool test_15() {
  Tensor<3,4,int> A, B;
  auto ait=A.begin();
  auto bit=B.begin();
  for(; !ait.end(); ait.next(), bit.next()) {
    *ait = 1;
    *bit = 1;
  }

  return contract<4>(A,B)==3*3*3*3;
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
  if(!test_9()) {
    retval |= 1<<8;
    std::cout << "Failed test_9" << "\n";
  }
  if(!test_10()) {
    retval |= 1<<9;
    std::cout << "Failed test_10" << "\n";
  }
  if(!test_11()) {
    retval |= 1<<10;
    std::cout << "Failed test_11" << "\n";
  }
  if(!test_12()) {
    retval |= 1<<11;
    std::cout << "Failed test_12" << "\n";
  }
  if(!test_13()) {
    retval |= 1<<12;
    std::cout << "Failed test_13" << "\n";
  }
  if(!test_14()) {
    retval |= 1<<13;
    std::cout << "Failed test_14" << "\n";
  }
  if(!test_15()) {
    retval |= 1<<14;
    std::cout << "Failed test_15" << "\n";
  }

  return retval;
}
