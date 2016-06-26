#ifndef __YAFEL_PCG_SOLVE_HPP
#define __YAFEL_PCG_SOLVE_HPP

/*
 * Implementation of the preconditioned conjugate-gradient iterative solver
 * for lienar systems with positive-definite, hermitian matrices.
 * The implementation follows the algorithm given in the wikipedia
 * article, which I believe follows Trefethen.
 */ 

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/access_sparse_matrix.hpp"
#include "lin_alg/operators.hpp"
#include "lin_alg/solver/iterative/Preconditioner.hpp"

YAFEL_NAMESPACE_OPEN

// not sure if I should use double or dataType for this.
// practically speaking, cg only makes sense for floating-point
// types, in which case double is fine.
constexpr double PCG_SOLVER_TOL = 1.0e-14;

template<typename MT, typename PT, typename dataType>
Vector<dataType> pcg_solve(const access_sparse_matrix<MT, dataType> &A,
                           const Vector<dataType> &rhs,
                           const Preconditioner<PT,dataType> &M,
                           dataType TOL = static_cast<dataType>(PCG_SOLVER_TOL))
{
    Vector<dataType> x0(rhs.size(), dataType(0));
    return pcg_solve(A,rhs,x0,TOL);
}


template<typename MT, typename PT, typename dataType>
Vector<dataType> pcg_solve(const access_sparse_matrix<MT, dataType> &A,
                           const Vector<dataType> &rhs,
                           const Preconditioner<PT,dataType> &M,
                           const Vector<dataType> &x0,
                           dataType TOL = static_cast<dataType>(PCG_SOLVER_TOL)) {
  
  
    Vector<dataType> x(x0);
    Vector<dataType> r = rhs - A*x0;
    Vector<dataType> z(r);
    M.solve(z);
  
    Vector<dataType> p(z);
  
    std::size_t k=0;
    std::size_t maxiter = x.size();
  
    dataType rTz_old = r.dot(z);
    dataType rTz_0 = rTz_old;
  
  
    if(rTz_old == dataType(0)) {
        return x; // <-- initial guess solved system exactly
    }
  
    while(k++ < maxiter) {

        auto Ap = A*p;
        dataType alpha = rTz_old/p.dot(Ap);
        x += p*alpha;
    
        r += -alpha*Ap;
        z = r;
        M.solve(z);
        
        dataType rTz_new = r.dot(z);

        //convergence check
        if(rTz_new/rTz_0 < TOL) {
            break;
        }

        dataType beta = rTz_new/rTz_old;
    
        rTz_old = rTz_new;
    
        p = z + beta*p;
    
    } // end while
  
    return x;
}






YAFEL_NAMESPACE_CLOSE


#endif
