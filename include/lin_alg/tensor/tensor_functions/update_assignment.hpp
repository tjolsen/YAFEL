//
// Created by tyler on 3/10/17.
//

#ifndef YAFEL_UPDATE_ASSIGNMENT_HPP_HPP
#define YAFEL_UPDATE_ASSIGNMENT_HPP_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"
#include "lin_alg/tensor/tensor_expression_types/BinaryOperations.hpp"

YAFEL_NAMESPACE_OPEN


template<typename T1, typename T2, int D, int R, typename dt1, typename dt2,
        bool b2, typename BinaryOperator>
inline T1 &update_assign(TensorExpression<T1, D, R, dt1, true> &&lhs,
                         const TensorExpression<T2, D, R, dt2, b2> &rhs, BinaryOperator)
{

    auto lit = lhs.begin();

    //automatically eval() the right-hand side for assignment.
    //going to have to be done anyway. Hopefully, the compiler
    //can see that the resulting values aren't going to be used anywhere,
    //so it'll just evaluate it on the fly
    for (auto r : rhs.eval()) {
        *lit = BinaryOperator::BinaryOp(*lit, r);
        ++lit;
    }

    return lhs.self();
}

template<typename T1, typename T2, int D, int R, typename dt1, typename dt2,
        bool b2, typename BinaryOperator>
inline T1 &update_assign(TensorExpression<T1, D, R, dt1, true> &lhs,
                         const TensorExpression<T2, D, R, dt2, b2> &rhs, BinaryOperator)
{
    auto lit = lhs.begin();

    //automatically eval() the right-hand side for assignment.
    //going to have to be done anyway. Hopefully, the compiler
    //can see that the resulting values aren't going to be used anywhere,
    //so it'll just evaluate it on the fly
    for (auto r : rhs.self().eval()) {
        *lit = BinaryOperator::BinaryOp(*lit, r);
        ++lit;
    }

    return lhs.self();
}


template<typename T1, int D, int R, typename dt1, typename dt2, typename BinaryOperator,
        typename=typename std::enable_if<ScalarTraits<dt2>::isYafelScalar()>::type>
inline T1 &update_assign(TensorExpression<T1, D, R, dt1, true> &lhs, dt2 rhs, BinaryOperator)
{

    for (auto &l : lhs) {
        l = BinaryOperator::BinaryOp(l, rhs);
    }

    return lhs.self();
}

template<typename T1, int D, int R, typename dt1, typename dt2, typename BinaryOperator,
        typename=typename std::enable_if<ScalarTraits<dt2>::isYafelScalar()>::type>
inline T1 &update_assign(TensorExpression<T1, D, R, dt1, true> &&lhs, dt2 rhs, BinaryOperator)
{

    for (auto &l : lhs) {
        l = BinaryOperator::BinaryOp(l, rhs);
    }

    return lhs.self();
}


//------------------------------------------------------------
// Update-assignment operators
//------------------------------------------------------------

template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b2>
auto operator+=(TensorExpression<T1, D, R, dt1, true> &lhs,
                const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return update_assign(lhs, rhs.self(), Addition<dt1, dt2>());
}

template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b2>
auto operator+=(TensorExpression<T1, D, R, dt1, true> &&lhs,
                const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return update_assign(std::forward<TensorExpression<T1, D, R, dt1, true>>(lhs), rhs.self(), Addition<dt1, dt2>());
}

template<typename T1, int D, int R, typename dt1, typename dt2,
        typename=typename std::enable_if<ScalarTraits<dt2>::isYafelScalar()>::type>
T1 &operator+=(TensorExpression<T1, D, R, dt1, true> &lhs, dt2 rhs)
{
    return update_assign(lhs, rhs, Addition<dt1, dt2>());
}

template<typename T1, int D, int R, typename dt1, typename dt2,
        typename=typename std::enable_if<ScalarTraits<dt2>::isYafelScalar()>::type>
T1 &operator+=(TensorExpression<T1, D, R, dt1, true> &&lhs, dt2 rhs)
{
    return update_assign(std::forward<TensorExpression<T1, D, R, dt1, true>>(lhs), rhs, Addition<dt1, dt2>());
}

//-----------------------------------------------------------------
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b2>
auto operator-=(TensorExpression<T1, D, R, dt1, true> &&lhs,
                const TensorExpression<T2, D, R, dt2, b2> &rhs)
{
    return update_assign(std::forward<TensorExpression<T1, D, R, dt1, true>>(lhs),
                         rhs, Subtraction<dt1, dt2>());
}

template<typename T1, int D, int R, typename dt1, typename dt2,
        typename=typename std::enable_if<ScalarTraits<dt2>::isYafelScalar()>::type>
T1 &operator-=(TensorExpression<T1, D, R, dt1, true> &&lhs, dt2 rhs)
{
    return update_assign(std::forward<TensorExpression<T1, D, R, dt1, true>>(lhs),
                         rhs, Subtraction<dt1, dt2>());
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UPDATE_ASSIGNMENT_HPP_HPP
