//
// Created by tyler on 2/26/17.
//

#ifndef YAFEL_TENSORCONTRACTION_HPP
#define YAFEL_TENSORCONTRACTION_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include "lin_alg/new_tensor/mp_utils/contraction_mp_utils.hpp"
#include <iostream>

YAFEL_NAMESPACE_OPEN

/**
 * \class TensorContraction
 * \brief Contractions of N indices of tensor expressions
 *
 * Note: Tensor contraction results in a tensor of rank R1 + R2 - 2*N
 *
 * @tparam TE1 TensorExpression type of lhs
 * @tparam TE2 TensorExpression type of rhs
 * @tparam D tensor dimension
 * @tparam R1 rank of lhs
 * @tparam R2 rank of rhs
 * @tparam N number of indices to contract
 * @tparam dt1 datatype of lhs
 * @tparam dt2 datatype of rhs
 * @tparam b1 "assignable" flag of lhs
 * @tparam b2 "assignable" flag of rhs
 */
template<typename TE1, typename TE2, int D, int R1, int R2, int N, typename dt1, typename dt2, bool b1, bool b2>
class TensorContraction : public TensorExpression<TensorContraction<TE1, TE2, D, R1, R2, N, dt1, dt2, b1, b2>, D,
        R1 + R2 - 2 * N, decltype(dt1(0) * dt2(0)), false>
{
public:
    using super = TensorExpression<TensorContraction<TE1, TE2, D, R1, R2, N, dt1, dt2, b1, b2>, D,
            R1 + R2 - 2 * N, decltype(dt1(0) * dt2(0)), false>;
    const TE1 &te_lhs;
    const TE2 &te_rhs;

    TensorContraction(const TensorExpression<TE1, D, R1, dt1, b1> &lhs,
                      const TensorExpression<TE2, D, R2, dt2, b2> &rhs)
            : te_lhs(lhs.self()), te_rhs(rhs.self()) {}


    decltype(dt1(0) * dt2(0)) linearIndexing(int idx) const noexcept
    {

        int lhs_offset = get_offset(typename super::stride_sequence(),
                                    get_first_N(sequence<R1 - N>(), typename TE1::super::stride_sequence()), idx);

        int rhs_offset = get_offset(typename super::stride_sequence(),
                                    get_last_N(sequence<(R1 == N) ? R2 - N : R2 - N + 1>(),
                                               typename TE2::super::stride_sequence()),
                                    idx % index_at(typename TE2::super::stride_sequence(), sequence<N - 1>()));

        return dot(
                make_const_slice(static_cast<typename TE1::super const &>(te_lhs),
                                 lhs_offset,
                                 get_last_N(sequence<N>(), typename TE1::super::stride_sequence())),
                make_const_slice(static_cast<typename TE2::super const &>(te_rhs),
                                 rhs_offset,
                                 get_first_N(sequence<N>(), typename TE2::super::stride_sequence()))
        );
    }


};

template<typename TE1, typename TE2, int D, int R1, int R2, int N, typename dt1, typename dt2, bool b1, bool b2>
inline auto contract(const TensorExpression<TE1, D, R1, dt1, b1> &lhs,
                     const TensorExpression<TE2, D, R2, dt2, b2> &rhs,
                     sequence<N>)
{
    return TensorContraction<TE1, TE2, D, R1, R2, N, dt1, dt2, b1, b2>(lhs, rhs);
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORCONTRACTION_HPP
