#ifndef __YAFEL_NULLELEMENT_HPP
#define __YAFEL_NULLELEMENT_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class NullElement : Element<NSD> {

public:
  using size_type = typename Element<NSD>::size_type;
  using coordinate_type = typename Element<NSD>::coordinate_type;


  NullElement(const DoFManager & dofm) : Element(dofm, ElementType::NULL_ELEMENT,
                                                 0, 0, 0, 0, 0)
  {}

  inline double shape_value_xi(size_type, const coordinate_type &) const {
    return 0;
  }

  inline double shape_grad_xi(size_type, size_type, const coordinate_type &) const {
    return 0;
  }

};


YAFEL_NAMESPACE_CLOSE

#endif
