//
// Created by tyler on 2/25/17.
//

#ifndef YAFEL_TENSOR_DOT_HPP
#define YAFEL_TENSOR_DOT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

/**
 * \brief Compute the scalar inner product of two equal rank/dimension tensors.
 *
 * Compute the full contraction of two tensors and return a scalar of type
 * decltype(dt1*dt2).
 * Both dt1 and dt2 must be initializable with zero for this to compile.
 * If you want a partial contraction of tensors, use the TensorContraction
 * expression instead.
 *
 * Note, this operation commutes. dot(A,B) == dot(B,A).
 * Do not confuse this with tensor contraction (of which matrix multiplication is a special case),
 * which does NOT commute, in general.
 *
 * @tparam T1 Tensor type of lhs
 * @tparam T2 Tensor type of rhs
 * @tparam D Tensor dimension
 * @tparam R Tensor rank
 * @tparam dt1 dataType of lhs
 * @tparam dt2 dataType of rhs
 * @tparam b1 "assignable" flag of lhs
 * @tparam b2 "assignable" flag of rhs
 * @param lhs Left-hand-side TensorExpression
 * @param rhs Right-hand-side TensorExpression
 * @return Full Contraction (scalar product) of lhs and rhs.
 */
template<typename T1, typename T2, int D, int R, typename dt1, typename dt2, bool b1, bool b2>
decltype(dt1()*dt2()) dot(const TensorExpression<T1,D,R,dt1,b1> &lhs,
         const TensorExpression<T2,D,R,dt2,b2> &rhs)
{

    auto lit=lhs.begin();
    auto rit=rhs.begin();
    decltype(dt1()*dt2()) retval(0);

    auto l_end = lhs.end();
    for(; lit != l_end ; ++lit, ++rit)
    {
        retval += (*lit)*(*rit);
    }

    return retval;
}


/**
 * \brief Compute the norm of a tensor
 *
 * Compute the L2 norm of a tensor by calling sqrt(dot(x,x))
 */
template<typename TE, int D, int R, typename dt, bool b>
auto norm(const TensorExpression<TE,D,R,dt,b> &x) {
    using std::sqrt;
    return sqrt(dot(x,x));
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOR_DOT_HPP
