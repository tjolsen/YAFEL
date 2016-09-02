#ifndef _YAFEL_LINTET_HPP
#define _YAFEL_LINTET_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class LinTet : public Element<NSD> {

public:
  using size_type = typename Element<NSD>::size_type;
  using coordinate_type = typename Element<NSD>::coordinate_type;

  LinTet(const DoFManager &dofm) : Element<NSD>(dofm, ElementType::LINEAR_TET,
                                           3, 1, 4*dofm.dof_per_node(), 10, 4)
  {
    double c = 1.0/4.0;
    this->xi_0.clear();
    this->quad_points.clear();
    this->quad_points.push_back(coordinate_type{c,c,c});
    this->quad_weights.resize(1, 1.0/6.0);

    this->xi_0.push_back(coordinate_type{0,0,0});
    this->xi_0.push_back(coordinate_type{1,0,0});
    this->xi_0.push_back(coordinate_type{0,1,0});
    this->xi_0.push_back(coordinate_type{0,0,1});
  }

  inline double shape_value_xi(size_type node, const coordinate_type &xi) const {
    switch(node) {
    case 0: return (1 - xi(0) - xi(1) - xi(2));
    case 1: return xi(0);
    case 2: return xi(1);
    case 3: return xi(2);
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
