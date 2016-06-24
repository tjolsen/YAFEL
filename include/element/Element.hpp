#ifndef _YAFEL_ELEMENT_HPP
#define _YAFEL_ELEMENT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/tensor_specializations.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/DoFManager.hpp"
#include "utils/ElementType.hpp"

#include <vector>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class Element {
  
public:
    using coordinate_type = Tensor<NSD,1>;
    using size_type = typename coordinate_type::size_type;


    size_type n_spaceDim; // number of spatial dimensions = dimension of coordinate_type
    size_type n_topoDim; // topological dimension of element (line = 1, quad/tri=2, hex/tet = 3)
    size_type n_quadPoints; // number of quadrature points. may be shifted into a QuadratureRule object in the future
    size_type dof_per_node; 
    size_type dof_per_el;
    int vtk_type;
    size_type el_num;
    size_type nodes_per_el;
  
    std::vector<double> quad_weights;
    std::vector<coordinate_type> quad_points;
    std::vector<coordinate_type> xi_0;
    std::vector<Tensor<NSD,2> > jacobians;
    std::vector<double> detJ;
    std::vector<Vector<double> > shape_vals; //vector of n_qp Vectors. holds values of shape funcs at qps
    std::vector<Matrix<double> > shape_grads; //vector of n_qp Matrix objects, hold grads wrt "x" of shape funcs at qps
    std::vector<Matrix<double> > shape_grads_0;//vector of n_qp Matrix objects, holds shape grads of canonical element at qps.
    std::vector<size_type> element_nodes;
    std::vector<coordinate_type> nodal_coords;
    std::vector<size_type> global_dofs;
    const DoFManager & DOFM;
    ElementType element_type;


    //  
    Element(const DoFManager &dofm, ElementType eltype, size_type ntd, 
            size_type nqp, size_type dofpe, int vtktype, size_type nodespe);
  
    // Virtual functions, specialized in child classes
    //virtual ~Element() {}

    virtual double shape_value_xi(size_type node, const coordinate_type &xi) const = 0;

    virtual double shape_grad_xi(size_type node, size_type component, const coordinate_type &xi) const = 0;

    //Functions in base class
    void calcJacobians(); // calcualte Jacobians at Gauss points and store in 
    void calcGrads(); // calculate shape function gradients (wrt spatial coords) and store in "grads[qpi](A, i)"
    void calcVals(); // calculate shape function values and store in "vals[qpi](A)"
    Tensor<NSD,2> calcJ_xi(const coordinate_type &xi) const;

    template<typename MTYPE, unsigned MNSD>
    void update_element(const GenericMesh<MTYPE,MNSD> & M, size_type elnum);

    inline double JxW(size_type qpi) const { return quad_weights[qpi]*detJ[qpi]; }
  
  
    //utility functions, might need to use in program so make public
    inline size_type getComp(size_type dof) const { return (dof % dof_per_node); }
    inline size_type getBase(size_type dof) const { return (dof/dof_per_node); }
  
    //return spatial interpolation at point xi
    coordinate_type xval(const coordinate_type & xi) const;

    
    // Function to initialize element. MUST be called before using element.
    // Cannot be placed into ctor due to virtual function calls.
    void init_element();
};


//================================================================
/*
 * Element implementation starts here
 */
//================================================================

//---------------------------------------------------------------------
template<unsigned NSD>
Element<NSD>::Element(const DoFManager &dofm, ElementType eltype, size_type ntd, size_type nqp, 
                      size_type dofpe, int vtktype, size_type nodespe) :
    n_spaceDim(NSD), 
    n_topoDim(ntd),
    n_quadPoints(nqp),
    dof_per_node(dofm.dof_per_node()),
    dof_per_el(dofpe), 
    vtk_type(vtktype), 
    nodes_per_el(nodespe),
    DOFM(dofm),
    element_type(eltype)
{

    //fill shape_grads vector of matrices
    shape_grads.clear();
    for(size_type qpi=0; qpi<n_quadPoints; ++qpi) {
        shape_grads.emplace_back(Matrix<double>(nodes_per_el, NSD, 0));
    }
}


//---------------------------------------------------------------------
template<unsigned NSD>
Tensor<NSD,2> Element<NSD>::calcJ_xi(const coordinate_type &xi) const {
  
    Tensor<NSD,2> ret;
  
    for(size_type i=0; i<NSD; ++i) {
        for(size_type j=0; j<NSD; ++j) {
            for(size_type A=0; A<nodes_per_el; ++A) {
                ret(i,j) += shape_grad_xi(A, j, xi) * nodal_coords[A](i);
            }
        }
    }
  
    return ret;
}

//---------------------------------------------------------------------
template<unsigned NSD>
void Element<NSD>::calcJacobians() {
    jacobians.clear();
    detJ.clear();
  
    for(size_type qpi=0; qpi<n_quadPoints; ++qpi) {

        Tensor<NSD,2> J;
        
        for(size_type i=0; i<NSD; ++i) {
            for(size_type j=0; j<NSD; ++j) {
                for(size_type A=0; A<nodes_per_el; ++A) {
                    J(i,j) += shape_grads_0[qpi](A, j) * nodal_coords[A](i);
                }
            }
        }

        jacobians.push_back(J);
    
        detJ.push_back(det(J));
    }
  
}


//---------------------------------------------------------------------
template<unsigned NSD>
void Element<NSD>::calcGrads() {
  
    for(size_type qpi=0; qpi<n_quadPoints; ++qpi) {
    
        auto Jinv = inv(jacobians[qpi]);
        auto & NG = shape_grads[qpi];

        for(size_type A=0; A<nodes_per_el; ++A) {
            for(size_type j=0; j<NSD; ++j) {
                NG(A,j) = 0;
                for(size_type k=0; k<NSD; ++k) {
                    NG(A,j) += shape_grads_0[qpi](A,k)*Jinv(k,j);
                }
            }
        }
    }
}


//---------------------------------------------------------------------
template<unsigned NSD>
void Element<NSD>::calcVals() {
    shape_vals.clear();
    for(size_type qpi=0; qpi<n_quadPoints; ++qpi) {
        auto qp = quad_points[qpi];
        Vector<double> V(nodes_per_el, 0.0);
        for(size_type A=0; A<nodes_per_el; ++A) {
            V(A) = shape_value_xi(A, qp);
        }
        shape_vals.push_back(V);
    }
}



//---------------------------------------------------------------------
template<unsigned NSD>
template<typename MTYPE, unsigned MNSD>
void Element<NSD>::update_element(const GenericMesh<MTYPE,MNSD> &M, size_type elnum) {

    // do not update a null element
    // ...you're going to have a bad time...
    if(element_type == ElementType::NULL_ELEMENT) {
        return;
    }
  
    size_type Nnodes = M.element(elnum).size();
  
    global_dofs.clear();
    nodal_coords.clear();
    element_nodes = M.element(elnum);
  
    for(size_type n=0; n<element_nodes.size(); ++n) {
    
        nodal_coords.push_back(M.node(element_nodes[n]));
        for(size_type j=0; j<dof_per_node; ++j) {
            size_type dofNum = DOFM.global_index(elnum, n, j);
            global_dofs.push_back(dofNum);
        }
    }
  
    calcJacobians();
    calcGrads();
    //calcVals();
}


template<unsigned NSD>
typename Element<NSD>::coordinate_type 
Element<NSD>::xval(const coordinate_type & xi) const {
  
    coordinate_type x;
    for(size_type A=0; A<nodes_per_el; ++A) {
        x += nodal_coords[A]*shape_value_xi(A,xi);
    }
  
    return x;
}



template<unsigned NSD>
void Element<NSD>::init_element() {

    //populate shape values
    shape_vals.clear();
    for(size_type qpi=0; qpi<n_quadPoints; ++qpi) {
        auto qp = quad_points[qpi];
        Vector<double> V(nodes_per_el, 0.0);
        for(size_type A=0; A<nodes_per_el; ++A) {
            V(A) = shape_value_xi(A, qp);
        }
        shape_vals.push_back(V);
    }
    
    //populate shape value gradients of parent element
    shape_grads_0.clear();
    for(size_type qpi=0; qpi<n_quadPoints; ++qpi) {
        shape_grads_0.emplace_back(Matrix<double>(nodes_per_el, NSD));

        for(size_type A=0; A<nodes_per_el; ++A) {
            for(size_type dim=0; dim<NSD; ++dim) {
                shape_grads_0[qpi](A,dim) = shape_grad_xi(A,dim,quad_points[qpi]);
            }
        }
    }

}

YAFEL_NAMESPACE_CLOSE

#endif //#ifndef __YAFEL_ELEMENT_HPP

