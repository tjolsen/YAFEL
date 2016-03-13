#include "LinearAdvection.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>

YAFEL_NAMESPACE_OPEN

//==============================================================
LinearAdvection::LinearAdvection(const AdvectionParameters &ap)
  : M(ap.mesh_dims, ap.dir_elems),
    AP(ap),
    dofm(M, 1),
    DGQ(ap.polyOrder, dofm, ap.Q2D, ap.Q1D),
    Me(DGQ.nodes_per_element, DGQ.nodes_per_element),
    Se(DGQ.nodes_per_element, DGQ.nodes_per_element),
    LU_Me(Me),
    vout(),
    vtkm(DGQ,M)
{
  //build mesh faces
  M.build_faces();
  
  //add mesh to vtk output object
  vout.addVTKObject(&vtkm);
  
}


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
    if(contract<1>(c,n) <= 0) {
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
  return c*((u_in + u_out)/2) - (1.0/2.0)*std::abs(contract<1>(c,n))*(u_out-u_in)*n;
  }
}

//==============================================================
void LinearAdvection::setup()
{
  
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

  std::vector<double> conserved_qty(Nt+1,0);
  conserved_qty[0] = integrate_field(u);
  
  for(std::size_t ti=1; ti<=Nt; ++ti) {
    
    Vector<double> k1 = AP.dt*residual(u);
    Vector<double> k2 = AP.dt*residual(u + 0.5*k1);
    Vector<double> k3 = AP.dt*residual(u + 0.5*k2);
    Vector<double> k4 = AP.dt*residual(u + k3);

    u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

    conserved_qty[ti] = integrate_field(u);
    write_output(u, ti);
  }

  //write conserved quantity out to a file
  std::string fname(AP.output_file_base);
  fname += "_conserved.csv";
  std::ofstream out(fname.c_str());
  for(std::size_t ti=0; ti<Nt; ++ti) {
    out << ti*AP.dt << ", " << conserved_qty[ti] << std::endl;
  }
  out.close();
}



//==============================================================
Vector<double> LinearAdvection::residual(const Vector<double> &u)
{
  
  std::size_t nnodes = DGQ.nodes_per_element;
  
  Vector<double> u_elem(nnodes);
  Vector<double> R(nnodes*M.n_elements());

  for(std::size_t elnum=0; elnum < M.n_elements(); ++elnum) {
    DGQ.update_nodes(M,elnum);
    DGQ.update_dofs(elnum);
    
    //passing reference to element to enable easier parallelization later.
    Vector<double> r_elem = element_residual(elnum, DGQ, u);
    
    //solve local residual (due to block-diag structure)
    Vector<double> MinvR_elem = LU_Me.linsolve(r_elem);

    //assemble local residual into global
    for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
      R(DGQ.global_dofs(A,0)) = MinvR_elem(A);
    }
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

  //stuff that won't change
  auto J = E.jacobians[0];//E.calc_J_xi(E.nodes_xi[A]);
  auto Jinv = inv(J);
  double detJ = det(J);

  
  //integrate faces
  std::size_t nfaces = 4;
  for(std::size_t f=0; f<nfaces; ++f) {

    auto n = E.mesh_face_normal(f);
    auto N = E.parent_face_normal(f);
    
    for(std::size_t fqpi=0; fqpi<E.Q1D.n_qp(); ++fqpi) {
      Face F = E.element_faces[f];
      //elnum of adjacent element
      std::size_t adj_elem = (elnum==F.inner) ? F.outer : F.inner;
      
      //face in adjacent element corresponding to interface
      std::size_t adj_elem_face = (elnum==F.inner) ? F.outer_face : F.inner_face;
        

      auto xi_qp = E.face_qp(f,fqpi);
      
      //compute u_in and u_out at quadrature point
      double u_in=0;
      double u_out=0;
      std::vector<double> vals(E.poly_order+1,0);
      for(std::size_t fA=0; fA<E.poly_order+1; ++fA) {
        std::size_t A_in = E.edge_nodes_ccw(fA,f);
        std::size_t A_out = E.edge_nodes_cw(fA,adj_elem_face);
        std::size_t A_out_global = adj_elem*nnodes + A_out;
        
        vals[fA] = E.shape_value_xi(A_in,xi_qp);
        u_in += u_elem(A_in)*vals[fA];
        u_out += u(A_out_global)*vals[fA];
      }
      for(std::size_t fA=0; fA<E.poly_order+1; ++fA) {
        
        auto flux =flux_function(u_in, u_out, AP.v_advection, n, F.boundary);
        r_elem(E.edge_nodes_ccw(fA,f)) -= vals[fA]*contract<1>(Jinv*flux,N)*detJ*E.face_weight(fqpi);
      }

      /*
      //compute flux*N at quadrature point
      for(std::size_t fA=0; fA<E.poly_order+1; ++fA) {
        //local node number
        std::size_t A = E.edge_nodes_ccw(fA, f);

        for(std::size_t fB=0; fB<E.poly_order+1; ++fB) {
          //local node number
          std::size_t B = E.edge_nodes_ccw(fB, f);
          
          //adjacent element local node number
          std::size_t adj_B = E.edge_nodes_cw(fB, adj_elem_face);
          
          //global dof of adjacent node
          std::size_t adj_global_dof = adj_elem*nnodes + adj_B;
          
          double u_in = u_elem(A);
          double u_out = u(adj_global_dof);
          
          auto flux = flux_function(u_in, u_out, AP.v_advection, n, F.boundary);
          //auto J = E.jacobians[0];//E.calc_J_xi(E.nodes_xi[A]);
          //auto Jinv = inv(J);
          //double detJ = det(J);
          
          r_elem(A) -= E.shape_value_xi(B,xi_qp)*
            (detJ*contract<1>(Jinv*flux,N)*E.shape_value_xi(A, xi_qp))*E.face_weight(fqpi);
        }
      }
      */
      
      /*      
      //distribute flux to nodes
      for(std::size_t fA=0; fA<E.poly_order+1; ++fA) {
        std::size_t A = E.edge_nodes_ccw(fA, f);
        r_elem(A) += qp_flux*E.shape_value_xi(A,xi_qp)*E.face_weight(fqpi);
      }
      */
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
void LinearAdvection::write_output(const Vector<double> &u, std::size_t ti)
{
  std::string fname(AP.output_file_base);
  char buf[10];
  sprintf(buf, "%06lu", ti);
  buf[6] = '\0';
  fname += "_";
  fname += buf;
  fname += ".vtu";

  std::cout << "writing ti=" << ti << " u(0)="<<u(0) <<std::endl;
  
  VTKScalarData vtku(u, VTKObject::VTKPOINTDATA, "u");
  vout.addVTKObject(&vtku);
  vout.write(fname);
  vout.clearData();
}


//==============================================================
double LinearAdvection::integrate_field(const Vector<double> &u)
{
  
  double integral = 0;

  for(std::size_t elnum=0; elnum<M.n_elements(); ++elnum) {
    //DGQ.update_nodes();
    DGQ.update_dofs(elnum);

    for(std::size_t qpi=0; qpi<DGQ.Q2D.n_qp(); ++qpi) {
      for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
	
	std::size_t GA = DGQ.global_dofs(A,0);
	
	integral += u(GA)*DGQ.shape_values[qpi](A)*DGQ.detJ[qpi]*DGQ.Q2D.weight(qpi);
      }
    }
  }
  return integral;
}


YAFEL_NAMESPACE_CLOSE
