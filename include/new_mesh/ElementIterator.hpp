//
// Created by tyler on 3/12/17.
//

#ifndef YAFEL_ELEMENTITERATOR_HPP
#define YAFEL_ELEMENTITERATOR_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

// Forward Declaration of Mesh
class Mesh;


/**
 * \class ElementIterator
 * \brief Iterator over nodes in an element
 */
class ElementIterator
{

public:

    ElementIterator(const Mesh &mesh, int element, int elemNode);


    /**
     * Get address of mesh. Useful for verifying identity of underlying mesh
     * between two iterators. Probably only used in debug mode.
     * @return
     */
    Mesh const* mesh_address() const noexcept
    {
        return &mesh_;
    }

    /**
     * Validate a pair of iterators. Used to ensure that:
     * 1) Iterators point to same mesh
     * 2) Iterators point to same element
     */
     bool valid_pair();


protected:
    const Mesh &mesh_;

};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ELEMENTITERATOR_HPP
