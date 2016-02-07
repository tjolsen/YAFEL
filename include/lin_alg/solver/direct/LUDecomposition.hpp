#ifndef _YAFEL_LUDECOMPOSITION_HPP
#define _YAFEL_LUDECOMPOSITION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"

YAFEL_NAMESPACE_OPEN

template<typename dataType>
class LUDecomposition {
  
  
public:
  using size_type = typename Matrix<dataType>::size_type;
  using value_type = typename Matrix<dataType>::value_type;
  using reference = typename Matrix<dataType>::reference;

  LUDecomposition(const Matrix<dataType> &A) : storage(A.rows(), A.rows(), dataType(0))
  {
    computeLU(A);
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
  
  void computeLU(const Matrix<dataType> &A) {
    size_type N = A.getRows();
    
    for(size_type i=0; i<N; ++i) {
      for(size_type j=0; j<N; ++j) {
        
        value_type sum = A(i,j);
        size_type lim = (i<j) ? i : j;
        
        for(size_type k=0; k<lim; ++k) {
          sum -=  L(i,k)*U(k,j);
        }
        
        if(i > j) {
          storage(i,j) = sum / U(j,j);
        }
        else {
          storage(i,j) = sum;
        }
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
  
    size_type N = storage.getRows();
    Vector<dataType> retVec(N, value_type(0));
    
    for(size_type i=N-1; i>=0; --i) {
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
