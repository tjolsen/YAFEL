#ifndef _YAFEL_DG_QUAD_HPP
#define _YAFEL_DG_QUAD_HPP

/*
 * Class for a quadrilateral element of polynomial order P
 * intended for use with discontinuous galerkin methods.
 *
 * Uses 2D lagrange polynomial interpolation, with basis nodes
 * laid out in a tensor-product ordering at points determined
 * by the Gauss-Lobatto quadrature rules. This minimizes the
 * Runge effect that plagues equal-spaced nodes at high order.
 *
 * Highly unstable, don't use unless you're me...
 */

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "lin_alg/tensor/tensor_specializations.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include "mesh/GenericMesh.hpp"
#include "mesh/Face.hpp"
#include "utils/ElementType.hpp"
#include "utils/DG_DoFManager.hpp"
#include "utils/DualNumber.hpp"
#include "utils/QuadratureRule.hpp"
#include "utils/GaussLobattoQuadrature.hpp"

#include <vector>
#include <exception>
#include <iostream>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD=2, typename dataType=double>
class DG_Quad {
  
public:
    using coordinate_type = Tensor<NSD,1,dataType>;
    using size_type = typename coordinate_type::size_type;
    using value_type = dataType;
  
    //number of solution variables at each node
    size_type dof_per_node;

    //polynomial interpolation order
    size_type poly_order;

    //nodes in element: = (poly_order + 1)^2
    size_type nodes_per_element;

    //Element Type
    ElementType element_type;
  
    //Reference to a DoFManager
    DoFManager & dofm;


    //locations of parent nodes in xi-space: [nNodes] x (NSD)
    std::vector<coordinate_type> nodes_xi;
  
    //locations of nodes in 1d. nodes_xi is formed by tensor product of this.
    std::vector<value_type> nodes_xi_1d;

    // (i,j) combinations of each node. Used for mapping linearized -> coordinate nodes.
    std::vector<Tensor<NSD,1,size_type> > nodes_ij;

    //locations of nodes in physical space: [nNodes] x (NSD)
    std::vector<coordinate_type> nodes_x;
  
    //shape function values at each node: [n_quadPoints] x [nNodes]
    std::vector<Vector<dataType> > shape_values;
  
    //shape function gradient at each node: [n_quadPoints] x (nNodes x NSD)
    std::vector<Matrix<dataType> > shape_gradients;
  
    //jacobians of mappings from parent to physical element: [nNodes] x (NSD x NSD)
    std::vector<Tensor<NSD,2,dataType> > jacobians;

    //determinants of jacobian mappings
    std::vector<dataType> detJ;
  
    //vector of global DoF numbers: (nNodes x dof_per_node)
    Matrix<size_type> global_dofs;
  
    //matrix holding node numbers of edges in a counter-clockwise direction: (poly_order+1 x 4)
    Matrix<size_type> edge_nodes_ccw;
  
    //matrix holding node numbers of edges in a clockwise direction: (poly_order+1 x 4)
    Matrix<size_type> edge_nodes_cw;

    //Mesh normals, assuming straight-edge elements
    std::vector<Tensor<NSD,1,dataType> > mesh_normals;

    //vector of element Faces (taken from mesh, after build_faces has been called)
    std::vector<Face> element_faces;
  
    // 2D quadrature rule for volume integrals
    QuadratureRule<NSD> &Q2D;

    // 1D quadrature rule for surface integrals
    QuadratureRule<NSD-1> &Q1D;
  
  
    /*
     * Constructor(s)
     */
    DG_Quad(size_type polyOrder, 
            DoFManager & _dofm, 
            QuadratureRule<NSD> &_q2d, 
            QuadratureRule<NSD-1> & _q1d);
  
  
    /*
     * Function that updates all element quantities
     */
    template<typename MT>
    void update_element(GenericMesh<MT,NSD> &M, size_type elnum);
  
  
    /*
     * Functions to fill members
     */
    template<typename MT>
    void update_nodes(const GenericMesh<MT,NSD> &M, size_type elnum);
  
    //template<typename MT>
    //void update_faces(const GenericMesh<MT,NSD> &M, size_type elnum);

    void update_dofs(size_type elnum);
    void calc_shape_values();
    void calc_jacobians();
    void calc_shape_gradients();


    /*
     * Utility functions that may be convenient for the user
     */
    //compute a shape function at any location in xi-space
    inline value_type shape_value_xi(size_type node, const coordinate_type & xi) const {
        return lagrange_poly_xi(xi, node, nodes_ij, nodes_xi_1d);
    }

    //compute the physical location of a node in xi-space
    coordinate_type xval(const coordinate_type &xi) const;
  
    //get physical normal vector to a face (assuming straight edges)
    inline Tensor<NSD,1,dataType> mesh_face_normal(size_type fnum) const {
        return mesh_normals[fnum];
    }

    //return number of dofs in element
    inline size_type dof_per_element() const { return dof_per_node*nodes_per_element; }

    //return base node of a dof
    inline size_type base(size_type dof) const {return dof/dof_per_node;}

    //return component number of a dof
    inline size_type component(size_type dof) const {return dof%dof_per_node;}

    //get xi-coordinate for face quadrature point
    inline Tensor<NSD,1,dataType> face_qp(size_type fnum, size_type fqpi) const {
        coordinate_type qp;
        switch(fnum) {
        case 0:
            return Tensor<NSD,1,dataType>{Q1D.qp(fqpi)(0), -1};
        case 1:
            return Tensor<NSD,1,dataType>{1, Q1D.qp(fqpi)(0)};
        case 2:
            return Tensor<NSD,1,dataType>{Q1D.qp(fqpi)(0), 1};
        case 3:
            return Tensor<NSD,1,dataType>{-1, Q1D.qp(fqpi)(0)};
        default:
            return Tensor<NSD,1,dataType>{0,0}; //<--garbage value not on an edge. don't wanna check for errors
        }
    }
  
    //get face quadrature weight
    inline dataType face_weight(size_type fqpi) const {
        return Q1D.weight(fqpi);
    }

    //get normal vector of parent (xi-space) element.
    //assumes that people aren't dumb, so doesn't look for out of bounds fnums
    //since the possibility of throwing an exception would slow it down probably
    inline Tensor<NSD,1,dataType> parent_face_normal(size_type fnum) const {
        switch(fnum) {
        case 0:
            return Tensor<NSD,1,dataType>{0, -1};
        case 1:
            return Tensor<NSD,1,dataType>{1, 0};
        case 2:
            return Tensor<NSD,1,dataType>{0, 1};
        case 3:
            return Tensor<NSD,1,dataType>{-1, 0};
        default:
            return Tensor<NSD,1,dataType>{0, 0}; //<--garbage value to catch errors
        }
    }
  
    /*
     * Utility functions that shouldn't be called by user,
     * but have been left public for the purposes of testing
     */
    /*
     *compute j-th lagrange polynomial at point xi using basis nodes
     *arranged via tensor product with, with xi_0 being the 1D points.
     *ij_pairs indicates the coordinate of the node:
   
     xi[node](dim) = xi_0(ij_pairs[node](dim))
   
     *Function is templated on T to allow for use with DualNumbers.
     */
    template<typename T>
    T lagrange_poly_xi(const Tensor<NSD,1,T> & xi, 
                       size_type node,
                       const std::vector<Tensor<NSD,1,size_type> > &ij_pairs,
                       const std::vector<value_type> & xi_0) const;

    template<typename T>
    T lagrange_poly_xi_1d(T xi, 
                          size_type j, 
                          const std::vector<value_type> & xi_0) const;
  
  
    /*
     * Compute jacobian transformation at a point
     */
    Tensor<NSD,2,dataType> calc_J_xi(const coordinate_type &xi) const;
  
}; //end class

//=====================================================================
/*
 * Implementation
 */
//=====================================================================
template<unsigned NSD,typename dataType>
DG_Quad<NSD,dataType>::DG_Quad(size_type polyOrder,
                               DoFManager & _dofm,
                               QuadratureRule<NSD> &_q2d,
                               QuadratureRule<NSD-1> &_q1d)
    : dof_per_node(_dofm.dof_per_node()),
      poly_order(polyOrder),
      nodes_per_element((poly_order+1)*(poly_order+1)),
      element_type(ElementType::DG_QUAD),
      dofm(_dofm),
      nodes_xi(),
      nodes_xi_1d(),
      nodes_ij(),
      nodes_x(),
      shape_values(),
      shape_gradients(),
      jacobians(),
      detJ(),
      global_dofs(nodes_per_element, dof_per_node),
      edge_nodes_ccw(poly_order+1, 4),
      edge_nodes_cw(poly_order+1, 4),
      mesh_normals(),
      element_faces(),
      Q2D(_q2d),
      Q1D(_q1d)
{
  
    if(polyOrder < 1) {
        throw std::invalid_argument("DG_Quad::Must use polynomial order of at least 1");
    }

    //lay out master element points according to gauss-lobatto quadrature points (see Hesthaven)
    GaussLobattoQuadrature<NSD> GLob(poly_order+1);
    nodes_xi = GLob.nodes;

    nodes_xi_1d = GLob.nodes_1d;
    nodes_ij = GLob.pairs;

    //fill shape_values
    calc_shape_values();


    //fill faces
    for(size_type i=0; i<polyOrder+1; ++i) {
        //face 0 (bottom face)
        edge_nodes_ccw(i,0) = i;
        edge_nodes_cw(i,0) = polyOrder - i;

        //face 1 (right face)
        edge_nodes_ccw(i,1) = (polyOrder+1)*i + polyOrder;
        edge_nodes_cw(i,1) = (polyOrder+1)*(polyOrder - i) + polyOrder;

        //face 2 (top face)
        edge_nodes_ccw(i,2) = (polyOrder+1)*polyOrder + (polyOrder-i);
        edge_nodes_cw(i,2) = (polyOrder+1)*polyOrder + i;
    
        //face 3 (left face)
        edge_nodes_ccw(i,3) = (polyOrder+1)*(polyOrder-i);
        edge_nodes_cw(i,3) = (polyOrder+1)*i;
    }

}

//------------------------------------------------------------------
template<unsigned NSD, typename dataType> template<typename MT>
void DG_Quad<NSD,dataType>::update_element(GenericMesh<MT,NSD> &M, size_type elnum)
{

    //faces must be built for this element. This is a one-time only cost.
    // Will not be called every time update_element is called
    if(!M.faces_built) {
        M.build_faces();
    }

    //update element node locations and face information
    update_nodes(M,elnum);

    //update face information
    //update_faces(M,elnum);

    //fill global dofs
    update_dofs(elnum);

    //fill jacobians at quadrature points
    calc_jacobians();

    //calculate gradients of shape functions
    calc_shape_gradients();
  
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType> template<typename MT>
void DG_Quad<NSD,dataType>::update_nodes(const GenericMesh<MT,NSD> &M, size_type elnum)
{

    //interpolate x-values of xi-nodes assuming (for now) straight-edged
    // bilinear interpolation from corners.
    std::vector<Tensor<NSD,1,size_type> > ij_corners{
        Tensor<NSD,1,size_type>{0,0},
            Tensor<NSD,1,size_type>{1,0},
                Tensor<NSD,1,size_type>{1,1},
                    Tensor<NSD,1,size_type>{0,1},
                        };
  
    std::vector<value_type> xi0_corners{-1, 1};

    std::vector<size_type> corner_mesh_nodes = M.element(elnum);
    std::vector<coordinate_type> x_corners;
    //restrict to grabbing first 4 coordinates, since we only want corners
    // and are assuming straight-edge interpolation
    for(size_type i=0; i<4; ++i) {
        x_corners.push_back(M.node(corner_mesh_nodes[i]));
    }

    //interpolate node locations
    nodes_x.clear();
    for(size_type A=0; A<nodes_per_element; ++A) {
        coordinate_type x;
        for(size_type B=0; B<x_corners.size(); ++B) {
            x += x_corners[B]*lagrange_poly_xi(nodes_xi[A], B, ij_corners, xi0_corners);
        }
        nodes_x.push_back(x);
    }

    //create mesh face normal vectors (cross prod of edge with 3-vector (0,0,1)
    mesh_normals.clear();
    for(size_type f=0; f<4; ++f) {
        coordinate_type dx = x_corners[(f+1)%4] - x_corners[f];
        coordinate_type n{dx(1), -dx(0)};
        n = n/std::sqrt(contract<1>(n,n));
        mesh_normals.push_back(n);
    }

    //fill element faces
    element_faces.clear();
    for(size_type f=0; f<4; ++f) {

        size_type fi = M.cell_faces[elnum][f];

        element_faces.push_back(M.mesh_faces[fi]);
    }

}

//--------------------------------------------------------
template<unsigned NSD, typename dataType>
void DG_Quad<NSD,dataType>::update_dofs(size_type elnum)
{
    for(size_type A=0; A<nodes_per_element; ++A) {
        for(size_type i=0; i<dof_per_node; ++i) {
            //hack: need to fix DG_DoFManager to play nice with this...
            //manually computing global dofs ASSUMING ONLY DG_QUADS IN MESH
            global_dofs(A,i) = (elnum*nodes_per_element + A)*dof_per_node + i;
        }
    }
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType>
void DG_Quad<NSD,dataType>::calc_shape_values()
{
    shape_values.clear();
    for(size_type qpi=0; qpi<Q2D.n_qp(); ++qpi) {
        shape_values.emplace_back(Vector<dataType>(nodes_per_element));

        for(size_type i=0; i<nodes_per_element; ++i) {
            shape_values[qpi](i) = lagrange_poly_xi(Q2D.qp(qpi), i, nodes_ij, nodes_xi_1d);
        }
    }
  
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType>
void DG_Quad<NSD,dataType>::calc_jacobians() {
    jacobians.clear();
  
    for(size_type qpi=0; qpi<Q2D.n_qp(); ++qpi) {
    
        auto J = calc_J_xi(Q2D.qp(qpi));
    
        jacobians.push_back(J);
        detJ.push_back(det(J));
    }//end qpi
  
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType>
void DG_Quad<NSD,dataType>::calc_shape_gradients() {
  
    shape_gradients.clear();
  
    for(size_type qpi=0; qpi<Q2D.n_qp(); ++qpi) {
    
        Tensor<NSD,2,dataType> Jinv = inv(jacobians[qpi]);
        Matrix<dataType> dNA_dxj(nodes_per_element, NSD);
        for(size_type A=0; A<nodes_per_element; ++A) {
            Tensor<NSD,1,DualNumber<dataType> > dual_xi;
            for(size_type i=0; i<NSD; ++i) {
                dual_xi(i) = make_dual(Q2D.qp(qpi)(i));
            }

            for(size_type j=0; j<NSD; ++j) {
                dual_xi(j).second = dataType(1);
                dataType dNAdxi_j = lagrange_poly_xi(dual_xi, A, nodes_ij, nodes_xi_1d).second;
        
                for(size_type i=0; i<NSD; ++i) {
                    dNA_dxj(A,i) += dNAdxi_j*Jinv(j,i);
                }
                dual_xi(j).second = dataType(0);
            }
      
        }//end A
        shape_gradients.push_back(dNA_dxj);
    }//end qpi
}


//--------------------------------------------------------
template<unsigned NSD, typename dataType>
typename DG_Quad<NSD,dataType>::coordinate_type
DG_Quad<NSD,dataType>::xval(const coordinate_type &xi) const
{
    coordinate_type x;
    for(size_type A=0; A<nodes_per_element; ++A) {
        x += shape_value_xi(A,xi)*nodes_x[A];
    }
  
    return x;
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType> template<typename T>
T DG_Quad<NSD,dataType>::lagrange_poly_xi(const Tensor<NSD,1,T> &xi, 
                                          size_type node,
                                          const std::vector<Tensor<NSD,1,size_type> > & ij_pairs,
                                          const std::vector<value_type> & xi_0) const
{
    T ret(1);
    for(size_type dim=0; dim<NSD; ++dim) {
        ret *= lagrange_poly_xi_1d(xi(dim), ij_pairs[node](dim), xi_0);
    }
    return ret;
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType> template<typename T>
T DG_Quad<NSD,dataType>::lagrange_poly_xi_1d(T xi, 
                                             size_type j,
                                             const std::vector<value_type> & xi_0) const
{
    T ret(1);
    for(size_type i=0; i<xi_0.size(); ++i) {
        if(i != j) {
            ret *= (xi - xi_0[i]) / (xi_0[j] - xi_0[i]);
        }
    }
    return ret;
}

//--------------------------------------------------------
template<unsigned NSD, typename dataType>
Tensor<NSD,2,dataType> DG_Quad<NSD,dataType>::calc_J_xi(const coordinate_type &xi) const
{

    Tensor<NSD,2,dataType> J;
    Tensor<NSD,1,DualNumber<dataType> > dual_xi;
    for(size_type i=0; i<NSD; ++i) {
        dual_xi(i) = make_dual(xi(i));
    }
    for(size_type A=0; A<nodes_per_element; ++A) {
        for(size_type j=0; j<NSD; ++j) {
            dual_xi(j).second = dataType(1);
      
            //shape func derivative wrt xi_j
            dataType dNAdxi_j = lagrange_poly_xi(dual_xi, A, nodes_ij, nodes_xi_1d).second;
      
            for(size_type i=0; i<NSD; ++i) {
                J(i,j) += dNAdxi_j*nodes_x[A](i);
            }
            dual_xi(j).second = 0;
        }
    }//end A

    return J;
}


YAFEL_NAMESPACE_CLOSE

#endif
