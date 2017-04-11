//
// Created by tyler on 2/17/17.
//

#ifndef YAFEL_ELEMENTTYPE_HPP
#define YAFEL_ELEMENTTYPE_HPP

#include "yafel_globals.hpp"
#include <tuple>

YAFEL_NAMESPACE_OPEN


/**
 * \class
 * \brief Enum to denote the topology of an element.
 *
 * TensorProduct --> {line/quad/hex}
 * Simplex --> {line/tri/tet}
 * None --> {} (used for template defaulting purposes)
 */
enum class ElementTopology : int8_t
{
    TensorProduct,
    Simplex,
    None
};

/**
 * \class ElementType
 *
 * Utility class to unambiguously define an element type.
 * Requires both the ElementClass (denoting topology),
 * a topological dimension "topoDim",
 * and a polynomial interpolation order "poly_order".
 *
 * An implementation of operator< is provided in order to
 * allow for use as a key in a std::map, if necessary.
 */
struct ElementType
{
    inline ElementType(ElementTopology ec=ElementTopology::None, int td=0, int po=0) : elementTopology(ec), topoDim(td), polyOrder(po)
    {}

    ElementTopology elementTopology;
    int8_t topoDim;
    int8_t polyOrder;

    inline bool operator<(const ElementType &rhs) const
    {
        return std::make_tuple(elementTopology, topoDim, polyOrder)
               < std::make_tuple(rhs.elementTopology, rhs.topoDim, rhs.polyOrder);
    }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ELEMENTTYPE_HPP
