#include "yafel_globals.hpp"
#include "element/DG_Quad.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "utils/GaussLegendreQuadrature.hpp"
#include "utils/GaussLobattoQuadrature.hpp"
#include "utils/DG_DoFManager.hpp"
#include "utils/DualNumber.hpp"
#include "lin_alg/tensor/Tensor.hpp"

#include <iostream>

using namespace yafel;


// create 1-element mesh on [0, 1]^2 and get the first detJ. (should equal 1/4)
// tons of stuff has to go right (and has already been debugged en route to this test)
// in order to get this correct.
bool test_1() {
  
    RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{1,1});
    M.build_faces();
    //use polynomials of order 15 along each edge, and 15-point 1d integration rule (super overkill!)
    //hopefully this high-order test exposes numerical stability issues, if present (which they are not)
    std::size_t N = 15;
    GaussLegendreQuadrature<2> Q2(N);
    GaussLegendreQuadrature<1> Q1(N);
    DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);
  
    DG_Quad<> DGQ(N, dofm, Q2, Q1);
    DGQ.update_element(M, 0);

    return std::abs(DGQ.detJ[0] - 1.0/4.0) < 1.0e-8;
}


/*
 * use shape_gradients to compute the gradient of a function (x^2 + y^2)
 * and compare with result obtained by using automatic differentiation
 */
template<typename T>
T func(T x, T y) {
    return x*x + y*y;
}

bool test_2() {

    RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{1,1});
    M.build_faces();
    std::size_t N = 2;
    GaussLegendreQuadrature<2> Q2(N);
    GaussLegendreQuadrature<1> Q1(N);
    DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

    DG_Quad<> DGQ(N, dofm, Q2, Q1);

    DGQ.update_element(M, 0);

    bool good = true;
    for(std::size_t qpi=0; qpi<DGQ.Q2D.n_qp(); ++qpi) {
    
        Tensor<2,1,double> xqp = DGQ.xval(DGQ.Q2D.qp(qpi));

        double f_x = func(DualNumber<double>(xqp(0),1), DualNumber<double>(xqp(1),0)).second;
        double f_y = func(DualNumber<double>(xqp(0),0), DualNumber<double>(xqp(1),1)).second;

        Tensor<2,1,double> gdual{f_x, f_y};
        Tensor<2,1,double> g;    
        for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
            auto x = DGQ.nodes_x[A];
            g(0) += func(x(0), x(1))*DGQ.shape_gradients[qpi](A,0);
            g(1) += func(x(0), x(1))*DGQ.shape_gradients[qpi](A,1);
       
        }

        double normdiff = contract<1>(g-gdual,g-gdual);
        good = good &&  normdiff < 1.0e-8;
    }
  
    return good;
}

//test faces built correctly. ccw and cw edges must be compatible
bool test_3() {
    RectilinearMesh<2> M(std::vector<double>{1,1}, std::vector<std::size_t>{2,2});
    M.build_faces();
    std::size_t N = 3;
    GaussLegendreQuadrature<2> Q2(N);
    GaussLegendreQuadrature<1> Q1(N);
    DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

    DG_Quad<> DGQ(N, dofm, Q2, Q1);
    DGQ.update_element(M, 0);

    bool good = true;
    for(std::size_t A=0; A<N+1; ++A) {
        for(std::size_t f=0; f<4; ++f) {
      
            good = good && (DGQ.edge_nodes_ccw(A,f)==DGQ.edge_nodes_cw(N-A,f));
        }
    }
  
    return true;
}


// integrate boundary of a rectangle to get perimeter
bool test_4() {

    double Lx = 1.5;
    double Ly = .25;
    RectilinearMesh<2> M(std::vector<double>{Lx, Ly}, std::vector<std::size_t>{1,1});
    M.build_faces();

    std::size_t N = 1;
    GaussLegendreQuadrature<2> Q2(N);
    GaussLegendreQuadrature<1> Q1(N);
    DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

    DG_Quad<> DGQ(N, dofm, Q2, Q1);
    DGQ.update_element(M, 0);
    double perimeter = 0;
    for(std::size_t f=0; f<4; ++f) {

        for(std::size_t eqpi=0; eqpi<DGQ.Q1D.n_qp(); ++eqpi) {      
            auto xi = DGQ.face_qp(f, eqpi);
            auto n = DGQ.mesh_face_normal(f);
            auto N = DGQ.parent_face_normal(f);
            auto J = DGQ.calc_J_xi(xi);
      
            auto Jinv = inv(J);
            auto detJ = det(J);
            perimeter += detJ*contract<1>(Jinv*n,N)*DGQ.face_weight(eqpi);
        }

    }//end eqpi
  
    return std::abs(perimeter - 2*(Lx + Ly)) < 1.0e-8;
}

//integrate perimeter of mesh using a multi-element mesh.
//will need to use the element_faces field of DG_Quad to do this.
bool test_5() {

    double Lx = 1.5;
    double Ly = .25;
    std::size_t Nx = 40;
    std::size_t Ny = 40;
    RectilinearMesh<2> M(std::vector<double>{Lx, Ly}, std::vector<std::size_t>{Nx,Ny});
    M.build_faces();
  
    std::size_t N = 2;
    GaussLegendreQuadrature<2> Q2(N);
    GaussLegendreQuadrature<1> Q1(N);
    DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

    DG_Quad<> DGQ(N, dofm, Q2, Q1);

    double perimeter = 0;

    for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {

        DGQ.update_element(M, elnum);
        for(std::size_t f=0; f<4; ++f) {
            if(!DGQ.element_faces[f].boundary) {
                continue; //this is an internal face. skip for this test.
            }
      
            for(std::size_t eqpi=0; eqpi<DGQ.Q1D.n_qp(); ++eqpi) {
                auto xi = DGQ.face_qp(f, eqpi);
                auto n = DGQ.mesh_face_normal(f);
                auto N = DGQ.parent_face_normal(f);
                auto J = DGQ.calc_J_xi(xi);
        
                auto Jinv = inv(J);
                auto detJ = det(J);
                perimeter += detJ*contract<1>(Jinv*n,N)*DGQ.face_weight(eqpi);
            }
      
        }//end eqpi
    }

    return std::abs(perimeter - 2*(Lx + Ly)) < 1.0e-8;
}


//integrate a constant "wind" blowing through mesh.
//should integrate to zero exactly, since no boundary
//conditions are changing it, and the constant wind
//is divergence-free
bool test_6() {

    //mesh parameters
    double Lx = 10;
    double Ly = 10;
    std::size_t Nx = 40;
    std::size_t Ny = 40;
    RectilinearMesh<2> M(std::vector<double>{Lx, Ly}, std::vector<std::size_t>{Nx,Ny});
    M.build_faces();

    //wind
    double vx = 1;
    double vy = 1;
    Tensor<2,1,double> v{vx,vy};
  
    std::size_t N = 2;
    GaussLegendreQuadrature<2> Q2(N);
    GaussLegendreQuadrature<1> Q1(N);
    DG_DoFManager<RectilinearMesh<2>,2> dofm(M, 1);

    DG_Quad<> DGQ(N, dofm, Q2, Q1);

    double integral = 0;

    for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {

        DGQ.update_element(M, elnum);
        for(std::size_t f=0; f<4; ++f) {
      
            for(std::size_t eqpi=0; eqpi<DGQ.Q1D.n_qp(); ++eqpi) {
                auto xi = DGQ.face_qp(f, eqpi);
                auto N = DGQ.parent_face_normal(f);
                auto J = DGQ.calc_J_xi(xi);
        
                auto Jinv = inv(J);
                auto detJ = det(J);
                integral += detJ*contract<1>(Jinv*v,N)*DGQ.face_weight(eqpi);
            }
      
        }//end eqpi
    }

    return std::abs(integral) < 1.0e-8;
}




int main() {

    int retval = 0;

    if(!test_1()) {
        std::cout << "Failed test_1()" << std::endl;
        retval |= 1<<0;
    }
    if(!test_2()) {
        std::cout << "Failed test_2()" << std::endl;
        retval |= 1<<1;
    }
    if(!test_3()) {
        std::cout << "Failed test_3()" << std::endl;
        retval |= 1<<2;
    }
    if(!test_4()) {
        std::cout << "Failed test_4()" << std::endl;
        retval |= 1<<3;
    }
    if(!test_5()) {
        std::cout << "Failed test_5()" << std::endl;
        retval |= 1<<4;
    }
    if(!test_6()) {
        std::cout << "Failed test_6()" << std::endl;
        retval |= 1<<5;
    }

    return retval;
}
