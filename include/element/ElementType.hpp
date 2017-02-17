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
enum class ElementClass {
    TensorProduct,
    Simplex,
    None
};


/**
 * \class ElementType
 *
 * Utility class to unambiguously define an element type.
 * Requires both the ElementClass (denoting topology)
 * and a polynomial interpolation order "poly_order".
 *
 * An implementation of operator< is provided in order to
 * allow for use as a key in a std::map, if necessary.
 */
struct ElementType {
    ElementClass element_class;
    int poly_order;

    inline bool operator<(const ElementType &rhs) {
        return std::make_tuple(element_class,poly_order)
               < std::make_tuple(rhs.element_class, rhs.poly_order);
    }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ELEMENTTYPE_HPP
