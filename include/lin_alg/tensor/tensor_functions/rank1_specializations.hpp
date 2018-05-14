//
// Created by tyler on 4/17/17.
//

#ifndef YAFEL_RANK1_SPECIALIZATIONS_HPP
#define YAFEL_RANK1_SPECIALIZATIONS_HPP


#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/tensor/tensor_functions/operators.hpp"
#include "utils/Range.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T1, typename T2, typename dT1, typename dT2, bool b1, bool b2>
auto cross(const TensorExpression<T1,3,1,dT1,b1> &lhs,
           const TensorExpression<T2,3,1,dT2,b2> &rhs)
{

    using retDT = decltype(std::declval<dT1>()*std::declval<dT2>());
    Tensor<3,1,retDT> retval;

    retval(0) = lhs(1)*rhs(2) - lhs(2)*rhs(1);
    retval(1) = lhs(2)*rhs(0) - lhs(0)*rhs(2);
    retval(2) = lhs(0)*rhs(1) - lhs(1)*rhs(0);
    return retval;
};

YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_RANK1_SPECIALIZATIONS_HPP
