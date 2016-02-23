#ifndef _YAFEL_LINTRI_HPP
#define _YAFEL_LINTRI_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class LinTri : public Element<NSD> {

public:
  using size_type = typename Element<NSD>::size_type;
  using coordinate_type = typename Element<NSD>::coordinate_type;

  LinTri(const DoFManager &dofm) : Element<NSD>(dofm, ElementType::LINEAR_TRI,
                                           2, 1, 3*dofm.getDofPerNode(), 5, 3)
  {
    double c = 1.0/3.0;
    this->xi_0.clear();
    this->quad_points.clear();
    this->quad_points.push_back(coordinate_type{c,c});
    this->quad_weights.resize(1, 1.0/2.0);

    this->xi_0.push_back(coordinate_type{0,0});
    this->xi_0.push_back(coordinate_type{1,0});
    this->xi_0.push_back(coordinate_type{0,1});
  }

  inline double shape_value_xi(size_type node, const coordinate_type &xi) const {
    switch(node) {
    case 0: return (1 - xi(0) - xi(1));
    case 1: return xi(0);
    case 2: return xi(1);
    }

    return 0;
  }
  inline double shape_grad_xi(size_type node, size_type component, const coordinate_type &) const {
    if(node==0) {
      return -1;
    }
    else {
      return (node-1)==component;
    }
  }
  
};

YAFEL_NAMESPACE_CLOSE

#endif
