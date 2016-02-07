#ifndef _YAFEL_ILUPRECONDITIONER_HPP
#define _YAFEL_ILUPRECONDITIONER_HPP

/*
 * Compue the ILU(0) factorization of a sparse matrix. Currently,
 * This is only implemented for the sparse_csr, but I will extend to
 * sparse_bcsr.
 *
 * The algorithm implemented herein may be found at:
 * http://www.cfd-online.com/Wiki/Incomplete_LU_factorization_-_ILU
 */


#include "yafel_globals.hpp"
#include "lin_alg/Preconditioner.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"

#include <cassert>
#include <iostream>

YAFEL_NAMESPACE_OPEN

template<typename dataType>
class ILUPreconditioner : public Preconditioner<ILUPreconditioner<dataType>, dataType> {
private:
  sparse_csr<dataType> ILU;
  
  using container_type = sparse_csr<dataType>;
  using value_type = typename sparse_csr<dataType>::value_type;
  using reference =  typename sparse_csr<dataType>::reference;
  using size_type = typename sparse_csr<dataType>::size_type;

public:
  // construct an incomplete LU factorization with no extra fill-in.
  // filling options could be added later, for S&G.
  // At first, no pivoting will be used, as this tends to increase the
  // bandwidth of resulting matrices arising from PDE numerical methods.
  ILUPreconditioner(const sparse_csr<dataType> &A) : ILU(A) {

#ifndef _OPTIMIZED
    assert(A.rows() == A.cols() && "ILUPreconditioner dimension mismatch");
#endif    

    for(size_type r=0; r<ILU.rows()-1; ++r) {
      bool flag = true;
      
      std::cout << r << std::endl;
      
      value_type d = value_type(1)/ILU(r, r, flag);
      
      for(size_type i=r+1; i<ILU.rows(); ++i) {

        reference A_ir = ILU(i,r,flag);
        
        if(!flag) continue; //skip if not in sparsity
        
        value_type e = d*A_ir;
        A_ir = e;
        
        size_type idxmin = ILU.row_ptr[i];
        size_type idxmax = ILU.row_ptr[i+1];
        for(size_type idx=idxmin; idx<idxmax; ++idx) {
          if(ILU.col_index[idx] > r) {
            idxmin = idx;
            break;
          }
        }

        for(size_type idx=idxmin; idx<idxmax; ++idx) {
          size_type j = ILU.col_index[idx];
          reference A_ij = ILU._data[idx];
          value_type A_rj = ILU(r,j,flag);
          
          if(flag) {
            A_ij -= e*A_rj;
          }

        } // end idx-loop
        
      } // end i-loop
      
    } // end r-loop
          
    
  } // end ctor

  Vector<dataType> MinvV(const Vector<dataType> &rhs) const;
  const container_type & getILU() {return ILU;}
};


YAFEL_NAMESPACE_CLOSE

#endif
