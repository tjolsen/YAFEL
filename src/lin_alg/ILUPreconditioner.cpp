#include "lin_alg/ILUPreconditioner.hpp"
#include <iostream>


YAFEL_NAMESPACE_OPEN

ILUPreconditioner::ILUPreconditioner(const sparse_csr &A) : ILU(A){

  unsigned N = ILU.getRows();
  
  // grab raw pointers to exploit sparsity
  unsigned *row_ptr = ILU.getRowPtr();
  unsigned *col_index = ILU.getColIndexPtr();
  double *data = ILU.getDataPtr();

  // Currently, no Pivoting is used, so be careful...
  
  bool flag = false;
  for(unsigned k = 0; k<N; ++k) {
    double dkk = 1.0/A(k,k,flag);
    if(!flag || dkk == 0.0) {
      std::cerr << "Warning, ILU failed" << std::endl;
    }
    for(unsigned i=k+1; i<N; ++i) {
      double aik = ILU(i,k,flag);
      if(!flag) {
	//(i,k) not in sparsity pattern, continue.
	continue;
      }
      
      double mik = aik*dkk;
      
      ILU.assign(i,k,mik);

      unsigned i_start=row_ptr[i];
      unsigned i_end = row_ptr[i+1];
      unsigned idx_start = i_start;
      for(unsigned idx = i_start; idx < i_end; ++idx) {
	unsigned j = col_index[idx];
	if(j > k) {
	  idx_start = idx;
	  break;
	}
      }
      
      for(unsigned idx = idx_start+1; idx<i_end; ++idx) {
	unsigned j = col_index[idx];
	data[idx] -= mik*ILU(k,j);
      }
      
    } // end i-loop
  } // end k-loop

}


Vector ILUPreconditioner::MinvV(const Vector &rhs) const {
  // setup
  unsigned N = ILU.getRows();
  Vector y(rhs);
  
  //grab ILU pointers in order to be able to exploit sparsity
  // grab raw pointers to exploit sparsity
  unsigned *row_ptr = ILU.getRowPtr();
  unsigned *col_index = ILU.getColIndexPtr();
  double *data = ILU.getDataPtr();

  // Two step solution to solving LU x = b
  // Step 1: L(Ux) = b, Ux = y;
  //   Forward substitute to obtain y=Ux
  // Step 2: Ux = y
  //   Backward substitution to obtain x = (LU)^{-1} * b
  
  // forward-substitute to obtain y=Ux = (L^1) * b
  for(unsigned i=0; i<N; ++i) {
    unsigned i_start = row_ptr[i];
    unsigned i_end = row_ptr[i+1];
    for(unsigned idx = i_start; idx<i_end; ++idx) {
      unsigned j = col_index[idx];
      if(j >= i) {
	//you have moved into the U half of the ILU
	break;
      }
      
      y(i) -= y(j)*data[idx];
      
    }//end idx-loop
  }//end i-loop


  // backward-substitute to obtain x = (LU)^{-1} * b
  for(unsigned i=0; i<N; ++i) {
    unsigned row = N-i-1;
    unsigned i_start = row_ptr[row];
    unsigned i_end = row_ptr[row+1];
    
    double diagval = 1;
    double diagidx = i_start;
    for(unsigned idx = i_start; idx<i_end; ++idx) {
      unsigned j = col_index[idx];
      
      if(j < row) {
	//still in L part of ILU
	continue;
      }
      else if(row == j) {
	//save info about the diag val, loc.
	diagval = data[idx];
	diagidx = idx;
	break;
      }
    }
    
    if(diagval == 0) {
      continue;
    }
    
    //idx loop broken apart to eliminate extra conditional checks
    //when we know how they will turn out
    for(unsigned idx=diagidx+1; idx<i_end; ++idx) {
      unsigned j = col_index[idx];
      y(row) -= y(j)*data[idx];
    }//end idx-loop

    y(row) /= diagval;
    
  }//end i-loop

  return y;
}


YAFEL_NAMESPACE_CLOSE
