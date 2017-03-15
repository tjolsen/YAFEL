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
enum class ElementClass : int
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
    inline ElementType(ElementClass ec, int td, int po) : elementClass(ec), topoDim(td), polyOrder(po)
    {}

    ElementClass elementClass;
    int topoDim;
    int polyOrder;

    inline bool operator<(const ElementType &rhs)
    {
        return std::make_tuple(elementClass, topoDim, polyOrder)
               < std::make_tuple(rhs.elementClass, rhs.topoDim, rhs.polyOrder);
    }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ELEMENTTYPE_HPP
