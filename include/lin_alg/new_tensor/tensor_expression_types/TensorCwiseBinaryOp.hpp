//
// Created by tyler on 3/8/17.
//

#ifndef YAFEL_TENSORCWISEBINARYOP_HPP
#define YAFEL_TENSORCWISEBINARYOP_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2, template<typename, typename> class BinaryOpType>
class TensorCwiseBinaryOp
        : public TensorExpression<TensorCwiseBinaryOp<T1, T2, D, R, dt1, dt2, b1, b2, BinaryOpType>, D, R, typename BinaryOpType<dt1, dt2>::result_type, false>
{
public:
    using super = TensorExpression<TensorCwiseBinaryOp<T1, T2, D, R, dt1, dt2, b1, b2, BinaryOpType>,
            D, R, typename BinaryOpType<dt1, dt2>::result_type, false>;


    using result_type = typename BinaryOpType<dt1, dt2>::result_type;

    const T1 &lhs_;
    const T2 &rhs_;

    TensorCwiseBinaryOp(const TensorExpression<T1, D, R, dt1, b1> &lhs,
                        const TensorExpression<T2, D, R, dt2, b2> &rhs)
            : lhs_(lhs.self()), rhs_(rhs.self()) {}

    result_type linearIndexing(int idx) const noexcept
    {
        return BinaryOpType<dt1, dt2>::BinaryOp(lhs_.linearIndexing(idx), rhs_.linearIndexing(idx));
    }


};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSORCWISEBINARYOP_HPP
