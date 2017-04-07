//
// Created by tyler on 2/17/17.
//

#ifndef YAFEL_ELEMENT_HPP
#define YAFEL_ELEMENT_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "ElementType.hpp"
#include "new_mesh/Mesh.hpp"
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
    Element(ElementType= {ElementTopology::None, 0, 0}, int dofPerNode = 1);

    // Struct that holds element type
    ElementType elementType;

    // Mesh of the local "Master" element
    Mesh localMesh;

    //Quadrature rule (ie, points and weights)
    QuadratureRule quadratureRule;

    // Shape function values and gradients (in parameter space)
    std::vector<Eigen::VectorXd> shapeValues;
    std::vector<Eigen::MatrixXd> shapeGradXi;

    // update element values at a quadrature point
    template<int NSD>
    void update(int elnum, int qpi, const DoFManager &dofm);

    //Element data at a quadrature point
    Eigen::MatrixXd shapeGrad;
    double detJ;
    std::vector<int> globalDofs;


    // useful getters for an element
    inline int dofPerNode() const { return dof_per_node; }

    inline int getNode(int dof) const { return dof / dof_per_node; }

    inline int getComp(int dof) const { return dof % dof_per_node; }

    inline int nQP() const { return static_cast<int>(quadratureRule.weights.size()); }

private:
    void make_simplex();

    void make_tensorProduct();

    int dof_per_node;
};


template<int NSD>
void Element::update(int elnum, int qpi, const DoFManager &dofm)
{

    std::vector<int> nodes;
    dofm.getGlobalNodes(elnum, nodes);

    Tensor<NSD,2> Jacobian(0);

    if(NSD == elementType.topoDim) {
        for (auto A : IRange(0, static_cast<int>(nodes.size()))) {
            for (auto j : IRange(0, elementType.topoDim)) {
                for (auto i : IRange(0, NSD)) {
                    Jacobian(i, j) += dofm.dof_nodes[nodes[A]](i)*shapeGradXi[qpi](A,j);
                }
            }
        }
    }

    //Tensor<NSD,2> Jinv =
}


YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_ELEMENT_HPP
