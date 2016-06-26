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
#include "lin_alg/solver/iterative/Preconditioner.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"

#include <cassert>
#include <iostream>

YAFEL_NAMESPACE_OPEN

template<typename dataType>
class ILUPreconditioner : public Preconditioner<ILUPreconditioner<dataType>, dataType> {  

public:
  using container_type = sparse_csr<dataType>;
  using value_type = typename sparse_csr<dataType>::value_type;
  using reference =  typename sparse_csr<dataType>::reference;
  using size_type = typename sparse_csr<dataType>::size_type;



  // construct an incomplete LU factorization with no extra fill-in.
  // filling options could be added later, for S&G.
  // At first, no pivoting will be used, as this tends to increase the
  // bandwidth of resulting matrices arising from PDE numerical methods.
  ILUPreconditioner(const container_type &A);


  void solve(Vector<dataType> &rhs) const;
  const container_type & getILU() {return ILU;}


private:
  container_type ILU;

  //backward substitution
  void b_subst(Vector<dataType> &rhs) const;
  
  //forward substitution
  void f_subst(Vector<dataType> &rhs) const;
};

//=============================================================================
/*
 * Implementation
 */
//=============================================================================

template<typename dataType>
ILUPreconditioner<dataType>::ILUPreconditioner(const sparse_csr<dataType> &A) 
  : ILU(A) 
{
#ifndef _OPTIMIZED
    assert(A.rows() == A.cols() && "ILUPreconditioner dimension mismatch");
#endif    

    for(size_type r=0; r<ILU.rows()-1; ++r) {
      bool flag = true;
      
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

}// end ctor


//---------------------------------------------------------

template<typename dataType>
void ILUPreconditioner<dataType>::solve(Vector<dataType> &rhs) const {
  
  f_subst(rhs); // <---
  b_subst(rhs); // <------- These functions modify the ret vector in place to eliminate excess copying

}


//---------------------------------------------------------

template<typename dataType>
void ILUPreconditioner<dataType>::f_subst(Vector<dataType> &rhs) const {
  
    for(size_type row=0; row<ILU.rows(); ++row) {

        size_type idxmin = ILU.row_ptr[row];
        size_type idxmax = ILU.row_ptr[row+1];
        
        for(size_type i=idxmin; i<idxmax; ++i) {

            size_type col = ILU.col_index[i];

            if(col >= row) {
                break;
            }

            rhs(row) -= ILU._data[i]*rhs(col);
        }
    }
    
}


//---------------------------------------------------------

template<typename dataType>
void ILUPreconditioner<dataType>::b_subst(Vector<dataType> &rhs) const {
    
    
    for(size_type idx=0; idx<ILU.rows(); ++idx) {
        
        size_type row = ILU.rows() - idx - 1;
        
        size_type idxmin = ILU.row_ptr[row];
        size_type idxmax = ILU.row_ptr[row+1];

        size_type i_diag(0);
        for(size_type i=idxmin; i<idxmax; ++i) {
            if(ILU.col_index[i] == row) {
                i_diag = i;
                break;
            }
        }

        for(size_type i=i_diag+1; i<idxmax; ++i) {
            size_type col = ILU.col_index[i];

            rhs(row) -= ILU._data[i]*rhs(col);
        }
        rhs(row) /= ILU._data[i_diag];
    }
    
}


YAFEL_NAMESPACE_CLOSE

#endif
