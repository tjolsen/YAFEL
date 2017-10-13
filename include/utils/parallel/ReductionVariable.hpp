//
// Created by tyler on 9/29/17.
//

#ifndef YAFEL_REDUCTIONVARIABLE_HPP
#define YAFEL_REDUCTIONVARIABLE_HPP

#include "yafel_globals.hpp"
#include "yafel_config.hpp"
#include <type_traits>

YAFEL_NAMESPACE_OPEN

namespace detail {

template<typename T>
using NeedsPadding = typename std::enable_if<sizeof(T) % config::cacheline_bytes != 0>::type;


constexpr int required_padding(int S) {
    auto spill = S % config::cacheline_bytes;
    return (spill == 0) ? 0 : config::cacheline_bytes - spill;
}

}//end namespace detail


/**
 * \class ReductionVariable
 * \brief Pads a type to an integer multiple of a cache line size
 * in order to prevent false sharing when using an array of them
 * across multiple threads in a shared memory environment.
 */
template<typename T, int RP = detail::required_padding(sizeof(T))>
struct alignas(alignof(T)) ReductionVariable : public T {
public:
    template<typename ...Args>
    ReductionVariable(Args && ...args) : T(std::forward<Args>(args)...) {}

private:
    int8_t padding[detail::required_padding(sizeof(T))];
};

/**
 * Specialization of ReductionVariable for types
 * that do not require padding
 */
template<typename T>
struct alignas(alignof(T)) ReductionVariable<T,0> : public T {};

YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_REDUCTIONVARIABLE_HPP
