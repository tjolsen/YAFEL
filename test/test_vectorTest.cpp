#include "old_handmade_linalg/Vector.hpp"
#include "old_handmade_linalg/operators.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>

using namespace yafel;


int main() {
  
  std::size_t N = 10;
  double val = 1.0;
  
  Vector<double> a(N,val);
  
  //Test proper length
  assert(a.size()==N && "TEST: Vector<double> constructed to correct size");

  // Test value populated correctly, and [] retrieves value
  double sum = 0;
  for(std::size_t i=0; i<a.size(); ++i) {
    sum += a(i);
  }
  assert( std::abs(sum - N*val)<1.0e-14 && "TEST: sum of elements correct" );


  // Test copy constructor
  Vector<double> b(a);
  assert(b==a && "TEST: Vector equality using copy constructor");

  // Test addition and scalar multiplication
  Vector<> c = a + b;
  Vector<> d = 2*a;
  assert(c==d && "TEST: Vector doubling using addition and scalar multiplication");

  // Test subtraction
  Vector<> e(N, 0);
  Vector<> f = a-b;
  assert(e==f && "TEST: Vector subtraction");

  return 0;
}
