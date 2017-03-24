#ifndef _YAFEL_ELEMENT_HPP
#define _YAFEL_ELEMENT_HPP


#include "yafel_globals.hpp"

#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/tensor_specializations.hpp"

#include "old_element/Element.hpp"
#include "old_element/ElementTraits.hpp"

#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"
#include "utils/DualNumber.hpp"

#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \brief Construct Element object based on element type EL
 */
template<typename EL>
Element make_element(const DoFManager &dofm) {
    using coordinate_type = Tensor<3,1,double>;
    
    Element E;
    E.DOFM = dofm;
    E.xi_0 = EL::parent_nodes();
    EL::quadrature(E.quad_points, E.quad_weights);
    
    
    std::vector<Matrix<double>> shape_grads(E.quad_points
    
    
    
    //Fill Traits
    E.element_type = ElementTraits<EL>::element_type;
    E.vtk_type = ElementTraits<EL>::vtk_type;
    E.nodes_per_el = ElementTraits<EL>::nodes_per_el;
    E.topo_dim = ElementTraits<EL>::DIM;

    //Allocate Element Storage Containers
    E.jacobians = std::vector<Matrix<double>>(E.quad_points.size(), 
                                              Matrix<double>(ElementTraits<EL>::DIM));
    E.element_nodes.resize(E.nodes_per_el,0);

}



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
