//
// Created by tyler on 5/23/17.
//

#ifndef YAFEL_POISSON_HPP
#define YAFEL_POISSON_HPP

#include "yafel_globals.hpp"
#include <eigen3/Eigen/Dense>

/**
 * Element matrix and local solve computation for
 * an HDG implementation of the Poisson equation.
 */
template<int NSD>
struct Poisson : PDEBase<NSD>
{
    //Stabilization constant
    constexpr static double tau = 1.0;
    constexpr static double k_inverse = 1.0;

    template<typename T>
    void UUTangent_boundary(const Element &E, int qpi,
                            const std::vector<int> &local_nodes,
                            Eigen::DenseBase<T> &UU);

    template<typename T>
    void UQTangent_volume(const Element &E, int qpi, Eigen::DenseBase<T> &UQ);

    template<typename T>
    void UQTangent_boundary(const Element &E, int qpi,
                            const std::vector<int> &local_nodes,
                            Eigen::DenseBase<T> &UQ);

    template<typename T>
    void UUhatTangent(const Element &E, int qpi, Eigen::DenseBase<T> &UUhat);

    template<typename T>
    void QUTangent_volume(const Element &E, int qpi, Eigen::DenseBase<T> &QU);

    template<typename T>
    void QQTangent(const Element &E, int qpi, Eigen::DenseBase<T> &QQ);

    template<typename T>
    void QUhatTangent(const Element &E, int qpi,
                      const std::vector<int> &local_nodes,
                      const std::vector<int> &face_nodes,
                      Tensor<NSD, 1> normal,
                      Eigen::DenseBase<T> &QUhat);

    template<typename T>
    void UhatUTangent(const Element &E, int qpi,
                      const std::vector<int> local_nodes,
                      const std::vector<int> face_nodes,
                      Eigen::DenseBase<T> &UhatU);

    template<typename T>
    void UhatQTangent(const Element &E, int qpi,
                      const std::vector<int> local_nodes,
                      const std::vector<int> face_nodes,
                      Tensor<NSD, 1> normal,
                      Eigen::DenseBase<T> &UhatQ);

    template<typename T>
    void UhatUhatTangent(const Element &E, int qpi, Eigen::DenseBase<T> &UhatUhat);

};


//---------------------------------------------
// Impelementation
//---------------------------------------------

template<int NSD>
template<typename T>
void Poisson<NSD>::UUTangent_boundary(const Element &E, int qpi,
                                      const std::vector<int> &local_nodes,
                                      Eigen::DenseBase<T> &UU)
{
    //This function is called in a "face integral" loop, so needs local_nodes.
    // Must be called after Element::face_update<NSD>, which will be called in HDGAssembly
    for (int i = 0; i < local_nodes.size(); ++i) {
        int A = local_nodes[i];
        for (int j = 0; j < local_nodes.size(); ++j) {
            int B = local_nodes[j];
            UU(A, B) += E.jxw * tau * E.boundaryShapeValues[qpi](i) * E.boundaryShapeValues[qpi](j);
        }
    }

}


template<int NSD>
template<typename T>
void Poisson<NSD>::UQTangent_volume(const Element &E, int qpi, Eigen::DenseBase<T> &UQ)
{
    for (int Anode = 0; Anode < UQ.rows(); ++Anode) {
        for (int Bnode = 0; Bnode < UQ.cols(); ++Bnode) {
            for (int i = 0; i < NSD; ++i) {
                int B = Bnode * NSD + i;
                UQ(Anode, B) += E.jxw * E.shapeValues[qpi](Anode) * E.shapeGrad(Bnode, i);
            }
        }
    }
}

template<int NSD>
template<typename T>
void Poisson<NSD>::UQTangent_boundary(const Element &E, int qpi,
                                      const std::vector<int> &local_nodes,
                                      Eigen::DenseBase<T> &UQ)
{
    //This function is called in a "face integral" loop, so needs local_nodes.
    // Must be called after Element::face_update<NSD>, which will be called in HDGAssembly
    for (int II = 0; II < local_nodes.size(); ++II) {
        int Anode = local_nodes[II];
        for (int JJ = 0; JJ < local_nodes.size(); ++JJ) {
            int Bnode = local_nodes[JJ];
            for (int i = 0; i < NSD; ++i) {
                int B = Bnode * NSD + i;
                UQ(Anode, B) -= E.jxw * tau * E.boundaryShapeValues[qpi](II) * E.boundaryShapeValues[qpi](JJ);
            }
        }
    }
}


template<int NSD>
template<typename T>
void Poisson<NSD>::UUhatTangent(const Element &E, int qpi,
                                const std::vector<int> &local_nodes,
                                const std::vector<int> &face_nodes,
                                Eigen::DenseBase<T> &UUhat)
{
    for (int i = 0; i < local_nodes.size(); ++i) {
        int A = local_nodes[i];
        for (int j = 0; j < face_nodes.size(); ++j) {
            int B = face_nodes[j];
            UUhat(A, B) -= E.jxw * tau * E.boundaryShapeValues[qpi](i) * E.boundaryShapeValues[qpi](j);
        }
    }
}


template<int NSD>
template<typename T>
void Poisson<NSD>::QUTangent_volume(const Element &E, int qpi, Eigen::DenseBase<T> &QU)
{
    for (int Anode = 0; Anode < E.shapeGrad.rows(); ++Anode) {
        for (int Bnode = 0; Bnode < E.shapeValues[qpi].rows(); ++Bnode) {
            for (int i = 0; i < NSD; ++i) {
                int A = Anode * NSD + i;
                QU(A, Bnode) += E.jxw * E.shapeGrade(Anode, i) * E.shapeValues[qpi](Bnode);

            }
        }
    }
}


template<int NSD>
template<typename T>
void Poisson<NSD>::QQTangent_volume(const Element &E, int qpi, Eigen::DenseBase<T> &QQ)
{
    int nNodes = E.localMesh.nNodes();
    for (int Anode = 0; Anode < nNodes; ++Anode) {
        for (int Bnode = 0; Bnode < nNodes; ++Bnode) {
            for (int i = 0; i < NSD; ++i) {
                int A = Anode * NSD + i;
                int B = Bnode * NSD + i;
                QQ(A, B) += k_inverse * E.jxw * E.shapeValues[qpi](Anode) * E.shapeValues[qpi](Bnode);
            }
        }
    }
}


template<int NSD>
template<typename T>
void Poisson<NSD>::QUhatTangent(const Element &E, int qpi,
                                const std::vector<int> &local_nodes,
                                const std::vector<int> &face_nodes,
                                Tensor<NSD, 1> normal,
                                Eigen::DenseBase<T> &QUhat)
{
    for (int I = 0; I < local_nodes.size(); ++I) {
        int Anode = local_nodes[I];
        for (int J = 0; J < face_nodes; ++J) {
            int Bnode = face_nodes[J];
            for (int i = 0; i < NSD; ++i) {

                int A = Anode * NSD + i;
                QUhat(A, Bnode) -= normal(i) * E.jxw * E.boundaryShapeValue[qpi](I) * E.boundaryShapeValue[qpi](J);

            }
        }
    }
}


template<int NSD>
template<typename T>
void Poisson<NSD>::UhatUTangent(const Element &E, int qpi,
                                const std::vector<int> local_nodes,
                                const std::vector<int> face_nodes,
                                Eigen::DenseBase<T> &UhatU)
{
    for (int I = 0; I < face_nodes.size(); ++I) {
        int Anode = face_nodes[I];
        for (int J = 0; J < local_nodes.size(); ++J) {
            int Bnode = local_nodes[J];
            UhatU(Anode, Bnode) -= tau * E.boundaryShapeValue[qpi](I) * E.boundaryShapeValue[qpi](J);
        }
    }
}


template<int NSD>
template<typename T>
void Poisson<NSD>::UhatQTangent(const Element &E, int qpi,
                                const std::vector<int> local_nodes,
                                const std::vector<int> face_nodes,
                                Tensor<NSD, 1> normal,
                                Eigen::DenseBase<T> &UhatQ)
{
    for (int I = 0; I < face_nodes.size(); ++I) {
        int Anode = face_nodes[I];
        for (int J = 0; J < local_nodes.size(); ++J) {
            int Bnode = local_nodes[J];
            for(int i=0; i<NSD; ++i) {
                int B = Bnode*NSD + i;
                UhatQ(Anode,B) = E.jxw*normal(i)*E.boundaryShapeValue[qpi](I)*E.boundaryShapeValue[qpi](J);
            }
        }
    }
}

#endif //YAFEL_POISSON_HPP
