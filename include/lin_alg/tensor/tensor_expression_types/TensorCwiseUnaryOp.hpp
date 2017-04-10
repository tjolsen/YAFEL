//
// Created by tyler on 3/9/17.
//

#ifndef YAFEL_TENSORCWISEUNARYOP_HPP
#define YAFEL_TENSORCWISEUNARYOP_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class TensorCwiseUnaryOp
 *
 * \brief Component-wise operation on a TensorExpression
 * @tparam TE TensorExpression type
 * @tparam D tensor Dimension
 * @tparam R tensor Rank
 * @tparam dt tensor dataType
 * @tparam b "assignable" flag
 */
template<typename TE, int D, int R, typename dt, bool b, template<typename> class UnaryOpType>
class TensorCwiseUnaryOp
        : public TensorExpression<TensorCwiseUnaryOp<TE, D, R, dt, b, UnaryOpType>, D, R, typename UnaryOpType<dt>::result_type, false>
{
public:
    using super = TensorExpression<TensorCwiseUnaryOp<TE, D, R, dt, b, UnaryOpType>, D, R, typename UnaryOpType<dt>::result_type, false>;
    using result_type = typename UnaryOpType<dt>::result_type;

    const TE & te_;

    TensorCwiseUnaryOp(const TensorExpression<TE,D,R,dt,b> &te)
            : te_(te.self())
    {}


    inline result_type linearIndexing(int idx) const noexcept
    {
        return UnaryOpType<dt>::UnaryOp(te_.linearIndexing(idx));
    }
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENWORCWISEUNARYOP_HPP
