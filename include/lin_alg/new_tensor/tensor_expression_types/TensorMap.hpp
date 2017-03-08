//
// Created by tyler on 3/8/17.
//

#ifndef YAFEL_TENSORMAP_HPP
#define YAFEL_TENSORMAP_HPP

#include "yafel_globals.hpp"
#include "lin_alg/new_tensor/TensorExpression.hpp"
#include "lin_alg/new_tensor/mp_utils/map_mp_utils.hpp"
#include <type_traits>

YAFEL_NAMESPACE_OPEN

/**
 * \class TensorMap
 * \brief Map an existing pre-allocated memory buffer as a TensorExpression.
 *
 * This class behaves the same as a Tensor, but rather than using
 * a stack-allocated internal buffer, it uses an externally allocated buffer.
 * Since it is constructed via pointer, it cannot check that enough memory
 * has been allocated.
 *
 * TensorMap respects const-ness of the underlying pointer type via
 * its metaprogramming PtrType and RefType constructs, which account
 * for whether this TensorMap is assignable.
 * Because of this complication, it is easier (and advisable) to construct
 * a TensorMap using the "make_TensorMap" function with an auto return type.
 * Eg:
 *
 * \code{.cpp}
 * double buffer[9];
 * auto myMap = make_TensorMap<3,2>(buffer);
 * \endcode
 *
 * This code block will map the "buffer" array as a dimension=3, rank=2
 * TensorMap, with the "double" type deduced by the compiler.
 * In this context, since "buffer" is not a pointer to const,
 * it can be assigned to via the usual Tensor machinery.
 *
 * \code{.cpp}
 * myMap(0,0) = 1;
 * \endcode
 *
 * @tparam D Tensor Dimension
 * @tparam R Tensor Rank/Order
 * @tparam dataType Data type of mapped buffer
 */
template<int D, int R, typename dataType, bool assignable>
class TensorMap : public TensorExpression<TensorMap<D, R, dataType, assignable>, D, R, dataType, assignable>
{
public:
    using PtrType = typename map_ptr<dataType, true, !assignable>::type;
    using RefType = typename map_ref<dataType, !assignable>::type;


    TensorMap(PtrType ptr) : dataPtr(ptr) {}

    template<typename=typename std::enable_if<assignable>::type>
    RefType linearIndexing(int idx) noexcept
    {
        return dataPtr[idx];
    }

    RefType linearIndexing(int idx) const noexcept
    {
        return dataPtr[idx];
    }

private:
    //Pointer to mapped memory
    PtrType dataPtr;

};

template<int D, int R, typename PTR>
auto make_TensorMap(PTR ptr)
{
    static_assert(std::is_pointer<PTR>::value, "Must construct TensorMap from a pointer");

    constexpr bool PtrToConst = std::is_const<PTR>::value;
    using value_type = typename std::remove_pointer<PTR>::type;

    return TensorMap<D, R, value_type, !PtrToConst>(ptr);
};

YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_TENSORMAP_HPP
