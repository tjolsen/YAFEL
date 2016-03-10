#include "LinearAdvection.hpp"

#include <cmath>
#include <iostream>

YAFEL_NAMESPACE_OPEN

//==============================================================
LinearAdvection::LinearAdvection(const AdvectionParameters &ap)
  : M(ap.mesh_dims, ap.dir_elems),
    AP(ap),
    dofm(M, 1),
    DGQ(ap.polyOrder, dofm, ap.Q2D, ap.Q1D),
    Me(DGQ.nodes_per_element, DGQ.nodes_per_element),
    Se(DGQ.nodes_per_element, DGQ.nodes_per_element),
    LU_Me(Me)
{}


//==============================================================
void LinearAdvection::run() {

  setup();

  //do time integration (handles output writing)
  RK4();
  
}

//==============================================================
Tensor<2,1,double> 
LinearAdvection::flux_function(double u_in, double u_out,
                               Tensor<2,1,double> c, Tensor<2,1,double> n,
                               bool is_boundary) const
{
  if(is_boundary) {
    //boundary conditions handled here
    if(contract<1>(c,n) < 0) {
      //inflow
      return Tensor<2,1,double>{0,0};
    }
    else {
      //outflow
      return u_in*c;
    }
  }
  else {
  //Upwind flux
  return c*((u_in + u_out)/2) + (1.0/2.0)*std::abs(contract<1>(c,n))*(u_out-u_in)*n;
  }
}

//==============================================================
void LinearAdvection::setup()
{
  
  M.build_faces();
  set_Me_Se();
  LU_Me.reinit(Me);

}

//==============================================================
void LinearAdvection::RK4()
{

  Vector<double> u = set_initial_condition();
  
  write_output(u, 0);
  
  //number of time steps
  std::size_t Nt = std::size_t(AP.Tfinal/AP.dt); //round down. no point being more precise now
  
  for(std::size_t ti=1; ti<=Nt; ++ti) {
    
    Vector<double> k1 = AP.dt*residual(u);
    Vector<double> k2 = AP.dt*residual(u + 0.5*k1);
    Vector<double> k3 = AP.dt*residual(u + 0.5*k2);
    Vector<double> k4 = AP.dt*residual(u + k3);

    u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    
    write_output(u, ti);
  }

}



//==============================================================
Vector<double> LinearAdvection::residual(const Vector<double> &u)
{
  
  std::size_t nnodes = DGQ.nodes_per_element;
  
  Vector<double> u_elem(nnodes);
  Vector<double> R(nnodes*M.n_elements());

  for(std::size_t elnum=0; elnum < M.n_elements(); ++elnum) {
    DGQ.update_element(M,elnum); // <--currently wasteful, don't update gradients to get more speed    
    
    //passing reference to element to enable easier parallelization later.
    Vector<double> r_elem = element_residual(elnum, DGQ, u);
    
    //solve local residual (due to block-diag structure)
    Vector<double> MinvR_elem = LU_Me.linsolve(r_elem);
  }

  return R;
}

//==============================================================
Vector<double> LinearAdvection::element_residual(std::size_t elnum, 
                                                 const DG_Quad<2,double> &E, 
                                                 const Vector<double> &u)
{
  std::size_t nnodes = E.nodes_per_element;
  Vector<double> u_elem(nnodes);

  //populate local solution from global vector
  for(std::size_t i=0; i<nnodes; ++i) {
    u_elem(i) = u( E.global_dofs(i,0) );
  }
  
  //put convection contribution into local residual
  Vector<double> r_elem = Se*u_elem;
  
  //integrate faces
  std::size_t nfaces = 4;
  for(std::size_t f=0; f<nfaces; ++f) {

    auto n = E.mesh_face_normal(f);
    auto N = E.parent_face_normal(f);
    
    for(std::size_t fqpi=0; fqpi<E.Q1D.n_qp(); ++fqpi) {

      auto xi_qp = E.face_qp(f,fqpi);
      
      double qp_flux = 0;

      //compute flux*N at quadrature point
      for(std::size_t fA=0; fA<E.poly_order+1; ++fA) {
        //local node number
        std::size_t A = E.edge_nodes_ccw(fA, f);
        
        Face F = E.element_faces[f];
        
        //elnum of adjacent element
        std::size_t adj_elem = (elnum==F.inner) ? F.outer : F.inner;

        //face in adjacent element corresponding to interface
        std::size_t adj_elem_face = (elnum==F.inner) ? F.outer_face : F.inner_face;
        
        //adjacent element local node number
        std::size_t adj_A = E.edge_nodes_cw(fA, adj_elem_face);

        //global dof of adjacent node
        std::size_t adj_global_dof = adj_elem*nnodes + adj_A;

        double u_in = u_elem(A);
        double u_out = u(adj_global_dof);

        auto flux = flux_function(u_in, u_out, AP.v_advection, n, F.boundary);
        auto J = E.calc_J_xi(E.nodes_xi[A]);
        auto Jinv = inv(J);
        double detJ = det(J);
        
        qp_flux += detJ*contract<1>(Jinv*flux,N)*E.shape_value_xi(A, xi_qp);
      }
      
      //distribute flux to nodes
      for(std::size_t fA=0; fA<E.poly_order+1; ++fA) {
        std::size_t A = E.edge_nodes_ccw(fA, f);
        r_elem(A) += qp_flux*E.shape_value_xi(A,xi_qp)*E.face_weight(fqpi);
      }
      
    }//end fqpi-loop
    
  } //end f-loop

  return r_elem;
}


//==============================================================
void LinearAdvection::set_Me_Se()
{
  
  //use first element to compute
  DGQ.update_element(M,0);

  Me *= 0;
  Se *= 0;

  for(std::size_t qpi=0; qpi<DGQ.Q2D.n_qp(); ++qpi) {
    double jxw = DGQ.detJ[qpi]*DGQ.Q2D.weight(qpi);

    for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
      for(std::size_t B=0; B<DGQ.nodes_per_element; ++B) {
        Me(A,B) += DGQ.shape_values[qpi](A)*DGQ.shape_values[qpi](B)*jxw;

        for(std::size_t i=0; i<2; ++i) {
          Se(A,B) += DGQ.shape_gradients[qpi](A,i)*AP.v_advection(i)*DGQ.shape_values[qpi](B)*jxw;
        }
      }
    }
  }

}

//==============================================================
Vector<double> LinearAdvection::set_initial_condition()
{
  
  Vector<double> u0(M.n_elements()*DGQ.nodes_per_element);

  for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {

    DGQ.update_element(M,elnum);
    
    for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
      u0(DGQ.global_dofs(A,0)) = AP.u0_func(DGQ.nodes_x[A]);
    }
    
  }
  
  return u0;
}

//==============================================================
void LinearAdvection::write_output(const Vector<double> &u, double ti)
{
  
  std::cout << "writing ti=" << ti << "u(0)=" << u(0) <<std::endl;
  

}



YAFEL_NAMESPACE_CLOSE
