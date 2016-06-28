#ifndef _YAFEL_ELEMENT_HPP
#define _YAFEL_ELEMENT_HPP


#include "yafel_globals.hpp"

#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/tensor_specializations.hpp"

#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"
#include "element/ElementTraits.hpp"

YAFEL_NAMESPACE_OPEN

/**
 * \class ElementImpl
 * PIMPL implementation class for the Element interface
 *
 * 
 */
class ElementImpl {
    
public:


    /**
     * 
     */
    template<unsigned DIM, typename T>
    virtual T shape_value_xi(size_type node, Tensor<DIM,T> xi) const = 0;
    //virtual double shape_value_xi(size_type node, const coordinate_type &xi) const = 0;
};


YAFEL_NAMESPACE_CLOSE


#endif
