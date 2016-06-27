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

    //loop over pivots
    for(size_type r=0; r<ILU.rows()-1; ++r) {

      
        //get pivot index and value
        size_type idx_pivot = ILU.row_ptr[r];
        value_type pivot(0);
        for(size_type idx=ILU.row_ptr[r]; idx<ILU.row_ptr[r+1]; ++idx) {
            size_type col = ILU.col_index[idx];
            if(col == r) {
                idx_pivot = idx;
                pivot = ILU._data[idx];
                break;
            }
        }
      
        if(pivot == 0) {
            std::cerr << "\\[\\e[1;31m\\]Warning: ILUPreconditioner: Found zero-valued pivot element.\n"
                      << "Continuing without eliminating with this pivot. \\[\\e[m\\]" << std::endl;
            continue;
        }

        //eliminate below pivot
        for(size_type i=r+1; i<ILU.rows(); ++i) {

            size_type idx_ir = ILU.row_ptr[i];
            bool in_sparsity(false);
            value_type A_ir(0);

            //get starting index in this row
            for(size_type idx = ILU.row_ptr[i]; idx<ILU.row_ptr[i+1]; ++idx) {
                size_type col = ILU.col_index[idx];
                if(col == r) {
                    idx_ir = idx;
                    in_sparsity = true;
                    A_ir = ILU._data[idx];
                }
                if(col > r) {
                    break;
                }
            }

            if(!in_sparsity) {
                continue;
            }

            
            //eliminate common elements in rows "r" and "i"
            value_type m_ir = A_ir/pivot;
            size_type idx_i = idx_ir;
            size_type idx_r = idx_pivot;
            size_type idx_imax = ILU.row_ptr[i+1];
            size_type idx_rmax = ILU.row_ptr[r+1];
            size_type col_r = ILU.col_index[idx_r];
            size_type col_i = ILU.col_index[idx_i];

            while(idx_i<idx_imax && idx_r<idx_rmax) {
                
                if(col_i == col_r) {
                    ILU._data[idx_i] -= m_ir*ILU._data[idx_r];
                    col_i = ILU.col_index[++idx_i];
                    col_r = ILU.col_index[++idx_r];
                }
                else if(col_i < col_r) {
                    col_i = ILU.col_index[++idx_i];
                }
                else if(col_i > col_r) {
                    col_r = ILU.col_index[++idx_r];
                }
                
            }
            ILU._data[idx_ir] = m_ir;

        }// end i-loop
      
    } // end r-loop (pivot loop)

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
