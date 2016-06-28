#ifndef _YAFEL_ELEMENTTRAITS_HPP
#define _YAFEL_ELEMENTTRAITS_HPP

#include "yafel_globals.hpp"
#include "utils/ElementType.hpp"
#include "utils/ElementVtkType.hpp"

YAFEL_NAMESPACE_OPEN

template<typename ElType>
struct ElementTraits {
    using size_type = std::size_t;

    static constexpr size_type DIM; ///< Topological Dimension of element (line=1, tri/quad=2, tet/hex=3)
    static constexpr size_type nodes_per_element; ///< Hopefully self-explanatory
    static constexpr ElementType eltype; ///< value from enum class ElementType (in include/utils/)
    static constexpr ElementVtkType vtktype; ///< value from enum class ElementVtkType (in include/utils/)
    
};


YAFEL_NAMESAPCE_CLOSE

#endif
