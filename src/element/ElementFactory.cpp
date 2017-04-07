//
// Created by tyler on 4/2/17.
//

#include "element/ElementFactory.hpp"


YAFEL_NAMESPACE_OPEN

ElementFactory::ElementFactory(int dof_per_node)
        : dof_per_node(dof_per_node), element_container()
{}

Element& ElementFactory::getElement(ElementType elementType)
{
    if(element_container.count(elementType) == 0) {
        auto eit = element_container.emplace(elementType, Element(elementType,dof_per_node));
        return (*(eit.first)).second;
    }
    else {
        return element_container[elementType];
    }
}


YAFEL_NAMESPACE_CLOSE