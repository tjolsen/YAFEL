#include "lin_alg/LUDecomposition.hpp"


YAFEL_NAMESPACE_OPEN

LUDecomposition::LUDecomposition(const FullMatrix &A) : storage(A.getRows())
{
  computeLU(A);
}

double LUDecomposition::L(unsigned i, unsigned j) const {
  return (i>j) ? storage(i,j) : ((i==j)?1.0 : 0.0);
}

double LUDecomposition::U(unsigned i, unsigned j) const {
  return (j>=i) ? storage(i,j) : 0.0;
}

void LUDecomposition::computeLU(const FullMatrix &A) {
  unsigned N = A.getRows();

  for(unsigned i=0; i<N; ++i) {
    for(unsigned j=0; j<N; ++j) {

      double sum = A(i,j);
      unsigned lim = (i<j) ? i : j;

      for(unsigned k=0; k<lim; ++k) {
	sum -=  L(i,k)*U(k,j);
      }

      if(i > j) {
	//L(i,j) = sum / U(j,j);
	storage(i,j) = sum / U(j,j);
      }
      else {
	//U(i,j) = sum;
	storage(i,j) = sum;
      }
    }
  }
  
  return;
}

// Solve linear system by forward-back substitution
Vector LUDecomposition::linsolve(const Vector &rhs) const {
  Vector y = f_subst(rhs);
  return b_subst(y);
}

// Function performing a forward-substitution solve with lower trianbular matrix L
Vector LUDecomposition::f_subst(const Vector & rhs) const {
  
  int N = storage.getRows();
  Vector retVec(N, 0.0);
  
  for(int i=0; i<N; ++i) {
    double val = rhs(i);
    for(int j=0; j<i; ++j) {
      val -= L(i,j)*retVec(j);
    }
    retVec(i) = val;
  }
  
  return retVec;
}

Vector LUDecomposition::b_subst(const Vector &rhs) const {
  
  int N = storage.getRows();
  Vector retVec(N, 0.0);
  
  for(int i=N-1; i>=0; --i) {
    double val = rhs(i);
    for(int j=i+1; j<N; ++j) {
      val -= U(i,j)*retVec(j);
    }
    retVec(i) = val/U(i,i);
  }
  
  return retVec;
}

YAFEL_NAMESPACE_CLOSE
