//
// Created by tyler on 2/25/17.
//

#ifndef YAFEL_TENSOR_DOT_HPP
#define YAFEL_TENSOR_DOT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
auto dot(const TensorExpression<T1,D,R,dt1,b1> &lhs,
         const TensorExpression<T2,D,R,dt2,b2> &rhs)
{

    auto lit=lhs.begin();
    auto rit=rhs.begin();
    decltype(dt1(0)*dt2(0)) retval(0);

    for(; lit != lhs.end(); ++lit, ++rit)
    {
        retval += (*lit)*(*rit);
    }

    return retval;
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOR_DOT_HPP
