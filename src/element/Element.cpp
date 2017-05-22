//
// Created by tyler on 3/14/17.
//

#include "element/Element.hpp"
#include "element/ShapeFunctionUtils.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

Element::Element(ElementType et, int dofpn)
        : elementType(et),
          localMesh(Mesh::DefinitionScheme::Explicit),
          dof_per_node(dofpn)
{

    switch (et.elementTopology) {
        case ElementTopology::TensorProduct:
            make_tensorProduct();
            break;
        case ElementTopology::Simplex:
            make_simplex();
            break;
        case ElementTopology::None:
            //null element. used for points/unsupported types
            break;
    }

    build_element_faces();
}


YAFEL_NAMESPACE_CLOSE