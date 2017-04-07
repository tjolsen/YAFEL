//
// Created by tyler on 4/2/17.
//

#ifndef YAFEL_ELEMENTFACTORY_HPP
#define YAFEL_ELEMENTFACTORY_HPP

#include "element/ElementType.hpp"
#include "element/Element.hpp"
#include <map>


YAFEL_NAMESPACE_OPEN

/**
 * \class ElementFactory
 */
class ElementFactory
{

public:
    ElementFactory(int dof_per_node = 1);

    Element &getElement(ElementType elementType);

    inline int getDofPerNode() const { return dof_per_node; }

private:
    int dof_per_node;
    std::map<ElementType, Element> element_container;

};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_ELEMENTFACTORY_HPP
