//
// Created by tyler on 2/25/17.
//

#ifndef YAFEL_TENSOR_DOT_HPP
#define YAFEL_TENSOR_DOT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/TensorExpression.hpp"
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
YAFEL_ALWAYS_INLINE decltype(dt1() * dt2()) dot(const TensorExpression <T1, D, R, dt1, b1> &lhs,
                            const TensorExpression <T2, D, R, dt2, b2> &rhs)
{

    decltype(dt1() * dt2()) s0(0), s1(0), s2(0), s3(0);
    constexpr auto N = T1::tensor_storage(R);
    constexpr int blk_size = 4;
    constexpr int Nblocks = N / blk_size;

    for (int i = 0; i < blk_size * Nblocks; i += blk_size) {
        s0 += lhs.linearIndexing(i) * rhs.linearIndexing(i);
        s1 += lhs.linearIndexing(i + 1) * rhs.linearIndexing(i + 1);
        s2 += lhs.linearIndexing(i + 2) * rhs.linearIndexing(i + 2);
        s3 += lhs.linearIndexing(i + 3) * rhs.linearIndexing(i + 3);
    }
    for (int i = blk_size * Nblocks; i < N; ++i) {
        s0 += lhs.linearIndexing(i) * rhs.linearIndexing(i);
    }

    return s0 + s1 + s2 + s3;

}


/**
 * \brief Compute the norm of a tensor
 *
 * Compute the L2 norm of a tensor by calling sqrt(dot(x,x))
 */
template<typename TE, int D, int R, typename dt, bool b>
auto norm(const TensorExpression <TE, D, R, dt, b> &x)
{
    using std::sqrt;
    return sqrt(dot(x, x));
};

/**
 * \brief Compute the squared norm of a tensor
 *
 * Compute the L2 squared norm of a tensor by calling dot(x,x)
 */
template<typename TE, int D, int R, typename dt, bool b>
auto norm2(const TensorExpression <TE, D, R, dt, b> &x)
{
    return dot(x, x);
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_TENSOR_DOT_HPP
