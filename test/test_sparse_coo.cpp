#include "yafel_globals.hpp"
#include "lin_alg/sparse_coo.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace yafel;


/*
 * test_1: 
 * - test bare construction from sequential additions,
 * - correct computation of rows(), cols(), and nnz()
 */
bool test_1() {
  
  std::size_t N = 10;
  sparse_coo<> C;
  for(std::size_t i=0; i<N; ++i) {
    C.add(i,i+1,1);
  }
  
  bool good_rows = (C.rows() == N);
  bool good_cols = (C.cols() == N+1);
  bool good_nnz = (C.nnz() == N);
  return good_rows && good_cols && good_nnz;
}


/*
 * test_2:
 * - test bare construction from sequential additions of identical elements
 * - correct computation of rows(), cols(), and nnz()
 */
bool test_2() {
  
  std::size_t N = 10;
  sparse_coo<> C;
  for(std::size_t i=0; i<N; ++i) {
    C.add(i,i,1);
  }
  for(std::size_t i=0; i<N; ++i) {
    C.add(i,i,1);
  }
  
  bool good_rows = (C.rows() == N);
  bool good_cols = (C.cols() == N);
  bool good_nnz = (C.nnz() == N);

  return good_rows && good_cols && good_nnz;
  
}


/*
 * test_3:
 * - test copy construction
 */
bool test_3() {

  std::size_t N = 10;
  sparse_coo<> C;
  for(std::size_t i=0; i<N; ++i) {
    C.add(i,i+1,1);
  }

  sparse_coo<> D(C);
  
  return C==D;
}

/*
 * test_4:
 * - test triplet construction
 */
bool test_4() {
  std::size_t N = 10;
  
  sparse_coo<> C;
  
  std::vector<typename sparse_coo<>::triplet> triplets;
  for(std::size_t i=0; i<N; ++i) {
    C.add(i,i+1,1);
    triplets.push_back( typename sparse_coo<>::triplet(i, i+1, 1) );
  }

  sparse_coo<> D(triplets);
  
  return C==D;
}


int main() {

  int retval = 0;

  if(!test_1()) {
    std::cout << "Failed test_1" << "\n";
    retval |= 1<<0;
  }

  if(!test_2()) {
    std::cout << "Failed test_2" << "\n";
    retval |= 1<<1;
  }

  if(!test_3()) {
    std::cout << "Failed test_3" << "\n";
    retval |= 1<<2;
  }
 
 if(!test_4()) {
    std::cout << "Failed test_4" << "\n";
    retval |= 1<<3;
  }
 
  
  return retval;
}
