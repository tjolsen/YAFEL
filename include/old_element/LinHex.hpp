#ifndef _YAFEL_LINHEX_HPP
#define _YAFEL_LINHEX_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"
#include "utils/DualNumber.hpp"
#include <vector>
#include <type_traits>
#include <cmath>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class LinHex : public Element<NSD> {

public:
    using size_type = typename Element<NSD>::size_type;
    using coordinate_type = typename Element<NSD>::coordinate_type;

    //template<typename = typename std::enable_if<NSD>=3>::type>
    LinHex(const DoFManager &dofm) : Element<NSD>(dofm, ElementType::LINEAR_HEX, 
                                                  3, 8, 8*dofm.dof_per_node(), 12, 8) {
        double a = 1.0/sqrt(3.0);
        this->quad_points.clear();
        this->quad_weights.resize(8,1.0);
        this->xi_0.clear();
    
        //set up parent element xi_0 vectors
        this->xi_0.emplace_back(coordinate_type{-1, -1, -1});
        this->xi_0.emplace_back(coordinate_type{ 1, -1, -1});
        this->xi_0.emplace_back(coordinate_type{ 1,  1, -1});
        this->xi_0.emplace_back(coordinate_type{-1,  1, -1});
        this->xi_0.emplace_back(coordinate_type{-1, -1,  1});
        this->xi_0.emplace_back(coordinate_type{ 1, -1,  1});
        this->xi_0.emplace_back(coordinate_type{ 1,  1,  1});
        this->xi_0.emplace_back(coordinate_type{-1,  1,  1});
    
        //assign gauss points and weights
        this->quad_points.emplace_back(coordinate_type{-a, -a, -a});
        this->quad_points.emplace_back(coordinate_type{ a, -a, -a});
        this->quad_points.emplace_back(coordinate_type{ a,  a, -a});
        this->quad_points.emplace_back(coordinate_type{-a,  a, -a});
        this->quad_points.emplace_back(coordinate_type{-a, -a,  a});
        this->quad_points.emplace_back(coordinate_type{ a, -a,  a});
        this->quad_points.emplace_back(coordinate_type{ a,  a,  a});
        this->quad_points.emplace_back(coordinate_type{-a,  a,  a});
    }

    
    template<typename T>
    T shape_value_xi_impl(size_type, Tensor<NSD,1,T>, std::false_type) const {
        return T(0);
    }
                 
    template<typename T>
    T shape_value_xi_impl(size_type node, Tensor<NSD,1,T> xi, std::true_type) const {
        return (this->xi_0[node](0)*xi(0) + 1)*(this->xi_0[node](1)*xi(1) + 1)*(this->xi_0[node](2)*xi(2) + 1) / 8.0;
    }

    double shape_value_xi(size_type node, const coordinate_type &xi) const override {
        return shape_value_xi_impl(node,xi, std::integral_constant<bool,NSD>=3>());
    }


    double shape_grad_xi(size_type node, size_type component, 
                         const coordinate_type &xi) const override {

        /*
          double val = 1.0;
          for(size_type i=0; i<this->n_topoDim; ++i) {
            if(component == i) {
                continue;
            }
            val *= (1.0/2.0) * ( this->xi_0[node](i)*xi(i) + 1 );
        }
        return val * (this->xi_0[node](component))/(2.0);
        */
        
        Tensor<NSD, 1, DualNumber<double>> xd;
        for(size_type i=0; i<NSD; ++i) {
            xd(i) = make_dual(xi(i));
        }
        xd(component).second = 1;
        auto val = shape_value_xi_impl(node,xd,std::integral_constant<bool,NSD>=3>());
        return val.second;
    }
  
};


YAFEL_NAMESPACE_CLOSE

#endif
