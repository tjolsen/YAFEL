//
// Created by tyler on 3/11/17.
//

#ifndef YAFEL_RANK2_SPECIALIZATIONS_HPP
#define YAFEL_RANK2_SPECIALIZATIONS_HPP

/**
 * \file
 * \brief Definition of functions/types useful for rank-2 tensors
 *
 * Contains:
 * - Trace
 * - operator* (for rank2*rank2 and rank2*rank1)
 * - Eye (D-dimensional rank-2 identity)
 * - Deviator
 * - Determinant (D=2, D=3)
 */

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/tensor/tensor_functions/operators.hpp"
#include "utils/Range.hpp"

YAFEL_NAMESPACE_OPEN

//=====================================================================
// Trace operator
//=====================================================================
template<typename TE, int D, typename dt, bool b>
dt trace(const TensorExpression<TE, D, 2, dt, b> &te)
{
    dt retval{0};
    for (auto i : IRange(0, D)) {
        retval += te(i, i);
    }
    return retval;
}


//=====================================================================
// Deviator
//=====================================================================
template<typename TE, int D, typename dt, bool b>
auto deviator(const TensorExpression<TE, D, 2, dt, b> &te)
{
    Tensor<D, 2, dt> ret(te);
    auto mean_trace = trace(te) / D;
    for (auto i : IRange(0, D)) {
        ret(i, i) -= mean_trace;
    }

    return ret;
}


//=====================================================================
// Symmetric part of tensor
//=====================================================================

template<typename TE, int D, typename dt, bool b>
auto sym(const TensorExpression<TE, D, 2, dt, b> &te)
{
    return (.5 * (te + te.template perm<1, 0>().eval())).eval(); // C++ treatment of dependent types is awkward...
}


//=====================================================================
// Skew symmetric part of tensor
//=====================================================================
template<typename TE, int D, typename dt, bool b>
Tensor<D, 2, dt> skw(const TensorExpression<TE, D, 2, dt, b> &te)
{
    return (.5 * (te - te.template perm<1, 0>())).eval(); // C++ treatment of dependent types is awkward...
}


//=====================================================================
// "multiplication" operator (contraction following convention)
// Only defined for rhs.rank()==1 or rhs.rank()==2
//=====================================================================
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2,
        typename=typename std::enable_if<R == 1 || R == 2>::type>
inline auto operator*(const TensorExpression<T1, D, 2, dt1, b1> &lhs,
                      const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return contract<1>(lhs, rhs);
}


//=====================================================================
// Identity Tensor
//=====================================================================
template<int D, typename dt=double>
inline Tensor<D, 2, dt> TensorEye()
{
    Tensor<D, 2, dt> ret;
    for (auto i : IRange(0, D))
        ret(i, i) = dt{1};

    return ret;
}


//===================================================================
// Tensor determinant
//===================================================================
template<typename dt=double>
inline dt determinant(const Tensor<1, 2, dt> &A)
{
    return A(0, 0);
}


template<typename dt=double>
inline dt determinant(const Tensor<2, 2, dt> &A)
{
    return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
}

template<typename dt=double>
inline dt determinant(const Tensor<3, 2, dt> &A)
{
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1))
           - A(0, 1) * (A(1, 0) * A(2, 2) - A(2, 0) * A(1, 2))
           + A(0, 2) * (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1));
}


//===================================================================
// Tensor inverse
//===================================================================
template<typename dt=double>
inline auto inverse(const Tensor<1, 2, dt> &A)
{
    Tensor<1, 2, dt> Ainv;
    Ainv(0,0) = 1.0/A(0,0);
    return Ainv;
};


template<typename dt=double>
inline Tensor<2, 2, dt> inverse(const Tensor<2, 2, dt> &A)
{
    dt detAinv = 1.0 / determinant(A);
    Tensor<2, 2, dt> Ainv;
    Ainv(0, 0) = A(1, 1) * detAinv;
    Ainv(1, 1) = A(0, 0) * detAinv;
    Ainv(0, 1) = -A(0, 1) * detAinv;
    Ainv(1, 0) = -A(1, 0) * detAinv;
    return Ainv;
};


template<typename dt=double>
inline Tensor<3, 2, dt> inverse(const Tensor<3, 2, dt> &A)
{
    auto detAinv = 1.0 / determinant(A);

    Tensor<3, 2, dt> Ainv;

    Ainv(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) * detAinv;
    Ainv(0, 1) = -(A(0, 1) * A(2, 2) - A(2, 1) * A(0, 2)) * detAinv;
    Ainv(0, 2) = (A(0, 1) * A(1, 2) - A(1, 1) * A(0, 2)) * detAinv;

    Ainv(1, 0) = -(A(1, 0) * A(2, 2) - A(2, 0) * A(1, 2)) * detAinv;
    Ainv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * detAinv;
    Ainv(1, 2) = -(A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) * detAinv;

    Ainv(2, 0) = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1)) * detAinv;
    Ainv(2, 1) = -(A(0, 0) * A(2, 1) - A(2, 0) * A(0, 1)) * detAinv;
    Ainv(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) * detAinv;

    return Ainv;
}

/**
 * Transpose a rank-2 tensor. In practice, the compiler
 * generates far superior code when "eval()" is called
 * after the perm<1,0>(). For this reason, it is put here
 * in the low-level code rather than leaving it to the user.
 * @return A^T
 */
template<int D, int R, typename dt>
Tensor<D,2,dt> transpose(Tensor<D,R,dt> &A)
{
    return (A.template perm<1, 0>()).eval();
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_RANK2_SPECIALIZATIONS_HPP
