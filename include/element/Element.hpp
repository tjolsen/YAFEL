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

#include <Eigen/Core>
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \class Element
 */
class Element
{

public:
    Element(ElementType et = {ElementTopology::None, 0, 0}, int dofPerNode = 1, int quadratureOrderMultiplier = 2);

    // Struct that holds element type
    ElementType elementType;

    // Mesh of the local "Master" element
    Mesh localMesh;
    std::vector<coordinate<>> boundaryNodes;

    //Quadrature rule (ie, points and weights)
    QuadratureRule quadratureRule;
    QuadratureRule boundaryQuadratureRule;

    //multiplier to select quadrature order for a given interpolation order.
    //Defaults to 2 because is useful for Galerkin FEM
    int quadratureOrderMultiplier;

    // Shape function values and gradients (in parameter space)
    std::vector<Eigen::VectorXd> shapeValues;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> shapeGradXi;

    std::vector<Eigen::VectorXd> boundaryShapeValues;
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> boundaryShapeGradXi;


    // Permutation arrays for face nodes
    //std::vector<std::vector<std::vector<std::vector<int>>>> face_perm;
    std::vector<std::vector<int>> face_nodes;

    //Element data at a quadrature point
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> shapeGrad;
    double jxw;
    std::vector<int> globalNodes;


    // update element values at a quadrature point
    template<int NSD>
    void update(int elnum, int qpi, const DoFManager &dofm);

    // update element values at a face quadrature
    // Returns the surface normal vector at that point
    //template<int NSD>
    //Tensor<NSD, 1> face_update(int elnum, int fqpi, const CellFace &F, const DoFManager &dofm);

    template<int NSD>
    Tensor<NSD, 1> face_update(int elnum, int fqpi, const std::vector<int> &local_fnodes, const DoFManager &dofm);

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

    void make_tensorProduct();

    void build_element_faces();

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
    dofm.getGlobalNodes(elnum, globalNodes);

    Tensor<NSD, 2> Jacobian(0);

    if (NSD == elementType.topoDim) {
        for (int A = 0; A < globalNodes.size(); ++A) {
            auto x = dofm.dof_nodes[globalNodes[A]];
            for(int i=0; i<NSD; ++i) {
                for (int j = 0; j < NSD; ++j) {
                    Jacobian(i, j) += x(i) * shapeGradXi[qpi](A, j);
                }
            }
        }

        //Add cases for (topoDim,NSD)==(2,3), (topoDim,NSD)==(1,2)
    }

    double detJ = determinant(Jacobian);
    jxw = detJ * quadratureRule.weights[qpi];
    Tensor<NSD,2> Jinv = inverse(Jacobian);
    Tensor<NSD, 2> JinvT = transpose(Jinv);

    if (shapeGrad.rows() != globalNodes.size() || shapeGrad.cols() != NSD) {
        shapeGrad.resize(globalNodes.size(), NSD);
    }

    //for (auto A : IRange(0, static_cast<int>(nodes.size()))) {
    for(int A=0; A<static_cast<int>(globalNodes.size()); ++A) {
        alignas(32) Tensor<NSD,1> tmpgrad = make_TensorMap<NSD, 1>(&shapeGradXi[qpi](A, 0));
        Tensor<NSD,1> tmpgrad2 = JinvT*tmpgrad;
        //for (auto d : IRange(0, NSD)) {
        for(int d=0; d<NSD; ++d) {
            //double s = dot(tmpgrad, JinvT(d,colon()));
            shapeGrad(A, d) = tmpgrad2(d);
        }
    }

}

template<int NSD>
Tensor<NSD, 1> Element::face_update(int elnum, int fqpi, const std::vector<int> &local_fnodes, const DoFManager &dofm)
{
    int fTopoDim = elementType.topoDim - 1;
    dofm.getGlobalNodes(elnum, globalNodes);

    // auto &local_fnodes = face_perm[flocal][frot][lr_flag];
    std::vector<int> fnodes;
    fnodes.reserve(local_fnodes.size());
    for (auto fn : local_fnodes) {
        fnodes.push_back(globalNodes[fn]);
    }

#ifndef NDEBUG
    if (fTopoDim == NSD) {
        throw std::runtime_error("Error: Using elements of higher topoDim than spatial dim.");
    }
#endif

    Tensor<NSD, 2> Jacobian;
    if (fTopoDim == NSD - 1) {
        //only going to use the first topoDim-1 columns of Jacobian
        for (int A = 0; A < fnodes.size(); ++A) {
            for (int i = 0; i < NSD; ++i) {
                for (int j = 0; j < NSD-1; ++j) {
                    Jacobian(i, j) += dofm.dof_nodes[fnodes[A]](i) * boundaryShapeGradXi[fqpi](A, j);
                }
            }
        }
    } else {
#ifndef NDEBUG
        throw std::runtime_error(
                std::string("Element::face_update: Invalid combination (NSD,fTopoDim) = (")
                + std::to_string(NSD)
                + ", "
                + std::to_string(fTopoDim)
                + ")"
        );
#endif
    }


    Tensor<NSD, 1> normal = getUnscaledNormal(Jacobian);
    double detJ = norm(normal);
    normal = normal / detJ;
    jxw = detJ * boundaryQuadratureRule.weights[fqpi];

    return normal;
}

YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_ELEMENT_HPP
