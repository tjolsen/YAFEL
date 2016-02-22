#ifndef _YAFEL_LINQUAD_HPP
#define _YAFEL_LINQUAD_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include <vector>
#include <type_traits>
#include <cmath>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class LinQuad : public Element<NSD> {

public:
  using size_type = typename Element<NSD>::size_type;
  using coordinate_type = typename Element<NSD>::coordinate_type;

  template<typename = typename std::enable_if<NSD>=2>::type>
  LinQuad(const DoFManager &dofm) : Element<NSD>(dofm, 2, 4, 4*dofm.getDofPerNode(), 9, 4) {
    double a = 1.0/sqrt(3.0);
    this->quad_points.clear();
    this->quad_weights.resize(4,1.0);
    this->xi_0.clear();
    
    //set up parent element xi_0 vectors
    this->xi_0.emplace_back(coordinate_type{-1, -1});
    this->xi_0.emplace_back(coordinate_type{ 1, -1});
    this->xi_0.emplace_back(coordinate_type{ 1,  1});
    this->xi_0.emplace_back(coordinate_type{-1,  1});
    
    //assign gauss points and weights
    this->quad_points.emplace_back(coordinate_type{-a, -a});
    this->quad_points.emplace_back(coordinate_type{ a, -a});
    this->quad_points.emplace_back(coordinate_type{ a,  a});
    this->quad_points.emplace_back(coordinate_type{-a,  a});
  }

  double shape_value_xi(size_type node, const coordinate_type &xi) const {
    double val = 1.0;
    for(unsigned i=0; i<this->n_topoDim; ++i) {
      val *= (1.0/2.0) * (this->xi_0[node](i)*xi(i) + 1);
    }
    
    return val;
  }

  double shape_grad_xi(size_type node, size_type component, const coordinate_type &xi) const {
    double val = 1.0;
    for(size_type i=0; i<this->n_topoDim; ++i) {
      if(component == i) {
        continue;
      }
      val *= (1.0/2.0) * ( this->xi_0[node](i)*xi(i) + 1 );
    }
    return val * (this->xi_0[node](component))/(2.0);
  }
  
};


YAFEL_NAMESPACE_CLOSE

#endif
