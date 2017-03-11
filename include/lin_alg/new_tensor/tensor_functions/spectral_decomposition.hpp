//
// Created by tyler on 3/11/17.
//

#ifndef YAFEL_SPECTRAL_DECOMPOSITION_HPP
#define YAFEL_SPECTRAL_DECOMPOSITION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_functions/rank2_specializations.hpp"
#include "utils/Range.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

/**
 *
 * \brief Compute the spectral decomposition of a real, symmetric
 * rank-2 tensor.
 *
 * This implements the Jacobi eigenvalue algorithm for real
 * symmetric tensors.
 * It uses a sequence of Givens rotations to reduce the size
 * of the off-diagonal elements.
 * The algorithm converges linearly, so long as the largest
 * off-diagonal is chosen as the pivot in each iteration.
 * It should be appropriate for D={2,3,4} tensors, but
 * larger dimension tensors should start to consider a
 * more sophisticated. (Eg: QR-based)
 *
 * The input to this function is assumed to be symmetric,
 * so No checking will be performed to verify this (since it
 * will almost exclusively be used by polar decompositions, which
 * can guarantee that it is called correctly.
 *
 * The original implementation of this was in C, written
 * for use in Abaqus (V)UMAT subroutines. It is generalized
 * and modernized herein to work natively with the TensorExpression
 * classes.
 *
 * \param A Input tensor. This will be modified and will end
 * up with eigenvalues on its diagonal
 *
 * \param Q Tensor with eigenvectors of the original A as columns
 */
template<int D, typename dt>
void spectral_decomposition(Tensor<D, 2, dt> &A, Tensor<D, 2, dt> &Q)
{
    using std::abs;

    // Convergence tolerance
    constexpr double TOL{1.0e-14};

    //Initialize eigenvectors to identity
    Q = TensorEye<D,dt>();
    Tensor<D,2,dt> A_old(A);
    Tensor<D,1,int> changed(1);



    while(true) { // <-- Isn't this exciting! Fortunately, this algorithm is "guaranteed" to converge...
        bool converged = true;
        for(auto c : changed)
        {
            converged = converged && (c==0);
        }
        if(converged) {
            break;
        }

        // Locate largest off-diagonal element
        int imax{0}, jmax{1};
        dt maxval{0};
        for (auto i : IRange(0, D)) {
            for (auto j : IRange(i + 1, D)) {
                if (std::abs(A(i, j)) > maxval) {
                    maxval = A(i, j);
                    imax = i;
                    jmax = j;
                }
            }
        }

        auto q = jacobi_rotation(A, imax, jmax);

        // apply rotation to A
        A = (q * A).eval() * q.template perm<1,0>(); // A <-- q*A*q'

        // update eigenvector rotation
        Q = (Q*q).eval();

        for(auto i : IRange(0,D)) {
            changed(i) = (abs(A(i,i) - A_old(i,i)) > TOL);
        }

        A_old = A;
    }
    Q = Q.template perm<1,0>().eval();
}


template<int dim, typename dt>
Tensor<dim, 2, dt> jacobi_rotation(const Tensor<dim, 2, dt> &D, int row, int col)
{
    // The "using" statements allow for potentially non-std types in dt (eg DualNumber),
    // for which sqrt() and abs() have been overloaded
    using std::sqrt;
    using std::abs;
    constexpr dt TOL = 1.0e-14; //Relative tolerance for something being "done"

    auto Q = TensorEye<dim, dt>();

    if (abs(D(row, col) / D(row, row)) < TOL) {
        //D(r,c) already eliminated. Return to avoid
        //a divide by zero
        return Q;
    }

    dt beta = (D(col, col) - D(row, row)) / (2 * D(row, col));
    dt t = -beta;
    auto sqrt_bb_p1 = sqrt(beta * beta + 1);
    if (beta > 0)
        t += sqrt_bb_p1;
    else {
        t -= sqrt_bb_p1;
    }

    dt c = dt{1} / sqrt(t * t + 1);
    dt s = t * c;

    Q(row, row) = c;
    Q(col, col) = c;
    Q(row, col) = -s;
    Q(col, row) = s;

    return Q;
}

YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_SPECTRAL_DECOMPOSITION_HPP
