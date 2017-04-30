//
// Created by tyler on 2/17/17.
//

#ifndef YAFEL_ELEMENT_HPP
#define YAFEL_ELEMENT_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "ElementType.hpp"
#include "mesh/Mesh.hpp"
#include "quadrature/QuadratureRule.hpp"
#include "utils/DoFManager.hpp"

#include <eigen3/Eigen/Core>
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \class Element
 */
class Element
{

public:
    Element(ElementType et = {ElementTopology::None, 0, 0}, int dofPerNode = 1);

    // Struct that holds element type
    ElementType elementType;

    // Mesh of the local "Master" element
    Mesh localMesh;
    std::vector<coordinate<>> boundaryNodes;

    //Quadrature rule (ie, points and weights)
    QuadratureRule quadratureRule;
    QuadratureRule boundaryQuadratureRule;

    // Shape function values and gradients (in parameter space)
    std::vector<Eigen::VectorXd> shapeValues;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> shapeGradXi;

    std::vector<Eigen::VectorXd> boundaryShapeValues;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> boundaryShapeGradXi;


    // Permutation arrays for face nodes
    std::vector<std::vector<std::vector<std::vector<int>>>> face_perm;

    //Element data at a quadrature point
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> shapeGrad;
    double jxw;
    std::vector<int> globalDofs;


    // update element values at a quadrature point
    template<int NSD>
    void update(int elnum, int qpi, const DoFManager &dofm);

    // update element values at a face quadrature
    // Returns the surface normal vector at that point
    template<int NSD>
    Tensor<NSD, 1> face_update(int elnum, int fqpi, const CellFace &F, const DoFManager &dofm);


    // useful getters for an element
    inline int dofPerNode() const { return dof_per_node; }

    inline int getNode(int dof) const { return dof / dof_per_node; }

    inline int getComp(int dof) const { return dof % dof_per_node; }

    //Get number of (volume) quadrature points
    inline int nQP() const { return static_cast<int>(quadratureRule.weights.size()); }

    //Get number of face/boundary quadrature points (per face)
    inline int nFQP() const { return static_cast<int>(boundaryQuadratureRule.weights.size()); }


private:
    void make_simplex();

    void build_tet_faces();

    void build_tri_faces();

    void build_quad_faces();

    void make_tensorProduct();

    inline Tensor<3, 1> getUnscaledNormal(Tensor<3, 2> Jacobian)
    {
        return cross(Jacobian(colon(), 0), Jacobian(colon(), 1));
    };

    inline Tensor<2, 1> getUnscaledNormal(Tensor<2, 2> Jacobian)
    {
        return Tensor<2, 1>{Jacobian(1, 0), -Jacobian(0, 0)};
    };

    int dof_per_node;
};


template<int NSD>
void Element::update(int elnum, int qpi, const DoFManager &dofm)
{

    std::vector<int> nodes;
    dofm.getGlobalNodes(elnum, nodes);

    Tensor<NSD, 2> Jacobian(0);

    if (NSD == elementType.topoDim) {
        for (auto A : IRange(0, static_cast<int>(nodes.size()))) {
            for (auto j : IRange(0, static_cast<int>(elementType.topoDim))) {
                for (auto i : IRange(0, NSD)) {
                    Jacobian(i, j) += dofm.dof_nodes[nodes[A]](i) * shapeGradXi[qpi](A, j);
                }
            }
        }

        //Add cases for (topoDim,NSD)==(2,3), (topoDim,NSD)==(1,2)
    }

    double detJ = determinant(Jacobian);
    jxw = detJ * quadratureRule.weights[qpi];
    Tensor<NSD, 2> Jinv = inverse(Jacobian);
    Tensor<NSD, 2> JinvT = Jinv.template perm<1, 0>();
    shapeGrad = Eigen::MatrixXd::Constant(nodes.size(), NSD, 0);

    for (auto A : IRange(0, static_cast<int>(nodes.size()))) {
        for (auto d : IRange(0, NSD)) {
            shapeGrad(A, d) = 0;
            double s = dot(make_TensorMap<NSD, 1>(&shapeGradXi[qpi](A, 0)), JinvT(d, colon()));
            shapeGrad(A, d) = s;
        }
    }


}

template<int NSD>
Tensor<NSD, 1> Element::face_update(int elnum, int fqpi, const CellFace &F, const DoFManager &dofm)
{
    int fTopoDim = elementType.topoDim - 1;
    std::vector<int> nodes;
    dofm.getGlobalNodes(elnum, nodes);

    int flocal{-1};
    int frot{-1};
    int lr_flag{-1};
    if (F.left == elnum) {
        flocal = F.left_flocal;
        frot = F.left_rot;
        lr_flag = 0;
    } else if (F.right == elnum) {
        flocal = F.right_flocal;
        frot = F.right_rot;
        lr_flag = 1;
    } else {
        throw std::runtime_error("Error: Non-corresponding elnum/CellFace in Element::face_update<NSD>");
    }

    auto &local_fnodes = face_perm[flocal][frot][lr_flag];
    std::vector<int> fnodes;
    fnodes.reserve(local_fnodes.size());
    for (auto fn : local_fnodes) {
        fnodes.push_back(nodes[fn]);
    }

    if (fTopoDim == NSD) {
        throw std::runtime_error("Error: Using elements of higher topoDim than spatial dim.");
    }

    Tensor<NSD, 2> Jacobian;
    if (fTopoDim == NSD - 1) {
        //only going to use the first topoDim-1 columns of Jacobian
        for (int A = 0; A < fnodes.size(); ++A) {
            for (int i = 0; i < NSD; ++i) {
                for (int j = 0; j < fTopoDim; ++j) {
                    Jacobian(i, j) += dofm.dof_nodes[fnodes[A]](i) * boundaryShapeGradXi[fqpi](A, j);
                }
            }
        }
    } else {
        throw std::runtime_error(
                std::string("Element::face_update: Invalid combination (NSD,fTopoDim) = (")
                + std::to_string(NSD)
                + ", "
                + std::to_string(fTopoDim)
                + ")"
        );
    }


    Tensor<NSD, 1> normal = getUnscaledNormal(Jacobian);
    double detJ = norm(normal);
    normal = normal / detJ;
    jxw = detJ * boundaryQuadratureRule.weights[fqpi];

    return normal;
}

template<int NSD>
Tensor<NSD, 1> getUnscaledNormal(Tensor<NSD, 2> Jacobian);

YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_ELEMENT_HPP
