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
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/new_tensor/tensor_functions/operators.hpp"
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
auto operator*(const TensorExpression<T1, D, 2, dt1, b1> &lhs,
               const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return contract<1>(lhs, rhs);
}


//=====================================================================
// Identity Tensor
//=====================================================================
template<int D, typename dt=double>
Tensor<3, 2, dt> TensorEye()
{
    Tensor<D, 2, dt> ret;
    for (auto i : IRange(0, D))
        ret(i, i) = dt{1};

    return ret;
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_RANK2_SPECIALIZATIONS_HPP