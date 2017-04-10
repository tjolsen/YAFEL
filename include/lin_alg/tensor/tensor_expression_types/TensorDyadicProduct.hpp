//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_TENSORDYADICPRODUCT_HPP
#define YAFEL_TENSORDYADICPRODUCT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T1, typename T2, int D, int R1, int R2,
        typename dt1, typename dt2, bool b1, bool b2>
class TensorDyadicProduct :
        public TensorExpression<TensorDyadicProduct<T1, T2, D, R1, R2, dt1, dt2, b1, b2>,
                D, R1 + R2, decltype(dt1() * dt2()), false>
{
public:
    using super = TensorExpression<TensorDyadicProduct<T1, T2, D, R1, R2, dt1, dt2, b1, b2>,
            D, R1 + R2, decltype(dt1() * dt2()), false>;
    using result_type = decltype(dt1()*dt2());

    const T1 &lhs_;
    const T2 &rhs_;

    TensorDyadicProduct(const TensorExpression<T1, D, R1, dt1, b1> &lhs,
                        const TensorExpression<T2, D, R2, dt2, b2> &rhs)
            : lhs_(lhs.self()), rhs_(rhs.self()) {}


    result_type linearIndexing(int idx) const noexcept
    {
        return lhs_.linearIndexing(lhs_index(idx, typename super::stride_sequence(), typename T1::super::stride_sequence()))
               *rhs_.linearIndexing(rhs_index(idx,typename super::stride_sequence()));
    }

private:
    template<int ...SS>
    inline int rhs_index(int idx, sequence<SS...>) const noexcept
    {
        return idx % index_at(sequence<SS...>(), sequence<R1 - 1>());
    }

    template<int ...SS>
    inline int lhs_index(int , sequence<SS...>, sequence<>) const noexcept
    {
        return 0;
    }

    template<int P, int ...PP, int S, int ...SS>
    inline int lhs_index(int idx, sequence<S, SS...>, sequence<P, PP...>) const noexcept
    {
        return P * (idx / S) + lhs_index(idx % S, sequence<SS...>(), sequence<PP...>());
    }

};

template<typename T1, typename T2, int D, int R1, int R2,
        typename dt1, typename dt2, bool b1, bool b2>
auto otimes(const TensorExpression<T1, D, R1, dt1, b1> &lhs,
            const TensorExpression<T2, D, R2, dt2, b2> &rhs)
{
    return TensorDyadicProduct<T1, T2, D, R1, R2, dt1, dt2, b1, b2>(lhs, rhs);
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORDYADICPRODUCT_HPP
