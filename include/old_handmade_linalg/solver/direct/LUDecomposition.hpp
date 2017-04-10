#ifndef _YAFEL_LUDECOMPOSITION_HPP
#define _YAFEL_LUDECOMPOSITION_HPP

#include "yafel_globals.hpp"
#include "old_handmade_linalg/sparse_csr.hpp"
#include "old_handmade_linalg/Matrix.hpp"
#include "old_handmade_linalg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

template<typename dataType>
class LUDecomposition {
  
  
public:
  using size_type = typename Matrix<dataType>::size_type;
  using value_type = typename Matrix<dataType>::value_type;
  using reference = typename Matrix<dataType>::reference;

  LUDecomposition(const Matrix<dataType> &A) : storage(A)
  {
    computeLU();
  }

  void reinit(const Matrix<dataType> &A) {
    storage = A;
    computeLU();
  }

  value_type L(size_type i, size_type j) const {
    return (i>j) ? storage(i,j) : ((i==j)?1.0 : 0.0);
  }

  value_type U(size_type i, size_type j) const {
    return (j>=i) ? storage(i,j) : 0.0;
  }

  Vector<dataType> linsolve(const Vector<dataType> & rhs) const {
    return b_subst(f_subst(rhs));
  }

private:
  Matrix<dataType> storage;
  
  void computeLU() {
    size_type N = storage.rows();
    
    for(size_type k=0; k<N; ++k) {
      for(size_type i=k+1; i<N; ++i) {
        value_type m = storage(i,k)/storage(k,k); //dies here if zero on diag
        for(size_type j=k; j<N; ++j) {
          storage(i,j) -= m*storage(k,j);
        }
        storage(i,k) = m;
      }
    }

    return;
  }
  
  Vector<dataType> f_subst(const Vector<dataType> &rhs) const {
    size_type N = storage.rows();
    Vector<dataType> retVec(N, value_type(0));
    
    for(size_type i=0; i<N; ++i) {
      value_type val = rhs(i);
      for(size_type j=0; j<i; ++j) {
        val -= L(i,j)*retVec(j);
      }
      retVec(i) = val;
    }
    
    return retVec;
  }
  
  Vector<dataType> b_subst(const Vector<dataType> &rhs) const {
  
    size_type N = storage.rows();
    Vector<dataType> retVec(N, value_type(0));

    for(size_type ii=1; ii<=N; ++ii) {
      size_type i = N - ii;
      value_type val = rhs(i);
      for(size_type j=i+1; j<N; ++j) {
        val -= U(i,j)*retVec(j);
      }
      retVec(i) = val/U(i,i);
    }
  
    return retVec;
  }

};

YAFEL_NAMESPACE_CLOSE

#endif
