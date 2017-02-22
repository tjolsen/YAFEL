//
// Created by tyler on 2/21/17.
//

#ifndef YAFEL_TENSORSUBTRACTION_HPP
#define YAFEL_TENSORSUBTRACTION_HPP

#include "lin_alg/new_tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename TE1, typename TE2, int D, int R, typename dataType>
class TensorSubtraction : public TensorExpression<TensorSubtraction<TE1, TE2, D, R, dataType>, D, R, dataType, false>
{

public:
    const TE1 &lhs;
    const TE2 &rhs;

    template<typename dt1, typename dt2, bool b1, bool b2>
    TensorSubtraction(const TensorExpression<TE1, D, R, dt1, b1> &l,
                      const TensorExpression<TE2, D, R, dt2, b2> &r)
            : lhs(l.self()), rhs(r.self()) {}


    inline dataType linearIndexing(int idx) const noexcept
    {
        return lhs.linearIndexing(idx) - rhs.linearIndexing(idx);
    }
};


template<typename TE1, typename TE2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
auto operator-(const TensorExpression<TE1, D, R, dt1, b1> &lhs,
               const TensorExpression<TE2, D, R, dt2, b2> &rhs) noexcept
{
    return TensorSubtraction<TE1, TE2, D, R, decltype(dt1(0) - dt2(0))>(lhs, rhs);
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORSUBTRACTION_HPP
