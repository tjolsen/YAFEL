#ifndef _YAFEL_LINLINE_HPP
#define _YAFEL_LINLINE_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "lin_alg/Vector.hpp"
#include "utils/DoFManager.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class LinLine : public Element<NSD> {

public:

    using size_type = typename Element<NSD>::size_type;
    using coordinate_type = typename Element<NSD>::coordinate_type;

    
    LinLine(const DoFManager &dofm)
	: Element<NSD>(dofm, ElementType::LINEAR_LINE, 1, 2, 2*dofm.dof_per_node(),3,2)
	{
	    //set up parent element xi_0 vectors
	    this->xi_0.clear();
	    //assign gauss points and weights
	    this->quad_points.clear();
	    this->quad_weights.clear();
	    
	    double a = 1.0/sqrt(3.0);
	    this->quad_points.push_back(coordinate_type{-a});
	    this->quad_points.push_back(coordinate_type{a});

	    this->quad_weights.resize(2,1.0);

	    this->xi_0.push_back(coordinate_type{-1});
	    this->xi_0.push_back(coordinate_type{1});
	}
    
    inline double shape_value_xi(size_type node, const coordinate_type &xi) const {

	if(node == 0)
	    return (1.0/2.0)*(1.0 - xi(0));
	else
	    return (1.0/2.0)*(xi(0) + 1);

    }

    inline double shape_grad_xi(size_type node, size_type , const coordinate_type &) const {
	if(node == 0)
	    return -1.0/2.0;
	else
	    return 1.0/2.0;
    }

    
};

YAFEL_NAMESPACE_CLOSE

#endif
