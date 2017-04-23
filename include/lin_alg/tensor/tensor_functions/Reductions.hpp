//
// Created by tyler on 4/23/17.
//

#ifndef YAFEL_REDUCTIONS_HPP_HPP
#define YAFEL_REDUCTIONS_HPP_HPP


#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"
#include "lin_alg/tensor/tensor_expression_types/BinaryOperations.hpp"


YAFEL_NAMESPACE_OPEN

/**
 * \brief perform a full reduction of a tensor expression using a supplied binary operation predicate
 *
 * @tparam T Tensor expression type
 * @tparam D tensor dimension
 * @tparam R tensor rank
 * @tparam dt tensor data type
 * @tparam b "assignable" bool for tensor expression
 * @tparam BinaryOpType Binary operation used for accumulating result
 * @param A Tensor expression on which to perform reduction
 * @return reduction value
 */
template<typename T, int D, int R, typename dt, bool b, template<typename, typename> class BinaryOpType>
auto full_reduction(const TensorExpression<T, D, R, dt, b> &A)
{
    auto accumulation = BinaryOpType<dt, dt>::identity_value();
    for (auto x : A) {
        accumulation = BinaryOpType<dt, dt>::BinaryOp(accumulation, x);
    }
    return accumulation;
};


template<typename T, int D, int R, typename dt, bool b>
auto sum(const TensorExpression<T, D, R, dt, b> &A)
{
    return full_reduction<T, D, R, dt, b, Addition>(A);
};

template<typename T, int D, int R, typename dt, bool b>
auto max(const TensorExpression<T, D, R, dt, b> &A)
{
    return full_reduction<T, D, R, dt, b, Max>(A);
};

template<typename T, int D, int R, typename dt, bool b>
auto min(const TensorExpression<T, D, R, dt, b> &A)
{
    return full_reduction<T, D, R, dt, b, Min>(A);
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_REDUCTIONS_HPP_HPP
