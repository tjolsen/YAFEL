#include "Poisson.hpp"
#include <iostream>

YAFEL_NAMESPACE_OPEN

Poisson::Poisson(const char *fname): 
    M(MeshReader::gmsh_read(std::string(fname))), DOFM(1),EF(M, DOFM),
    BC(M, DOFM, 1, 0, SpatialFunction<double>([](const Vector &x){return 1.0;})),
    BC2(M,DOFM, 2, 0, SpatialFunction<double>([](const Vector &x){return 9.0;}))
{}

void Poisson::setup() {
  
    Fsys = Vector(EF.get_n_dof(), 0.0);


    // Boundary condition setup moved to constructor for now  
    /*
      bcnodes.clear();
      bcvals.clear();

      //set up boundary conditions: zero on boundary (physical id = 1)
      unsigned DIRICHLET_BOUNDARY = 1;
  
      for(unsigned elnum=0; elnum<M.get_n_elems(); ++elnum) {
      if(M.el_tags[elnum][0] == DIRICHLET_BOUNDARY) {
      for(unsigned n=0; n<M.elements[elnum].size(); ++n) {
      bcnodes.push_back(M.elements[elnum][n]);
      bcvals.push_back(0.0);
      }
      }
      }
    */

    //set volumetric source term
    fvol = 1;

    std::cout << "setup done\n";
}


void Poisson::assemble() {
  
    //loop over all elements
    std::cout << "Assembling...\n";
    for(unsigned elnum=0; elnum < M.get_n_elems(); ++elnum) {
    
        Element *e = EF.getElement(elnum);
        if(e == NULL) {
            continue;
        }
        if(e->n_spaceDim != 2) {
            continue;
        }
        //this next block is specific to testing the "cracked" mesh cracked_simpleCracked.msh
        if(M.el_tags[elnum][0] == 11) {
            continue;
        }

        // update element shape function values, gradients, etc
        e->update_element(M, elnum);
        int Ndof = e->dof_per_el;
    
        FullMatrix Kel(Ndof, Ndof, 0);
        Vector Fel(Ndof, 0.0);
    
        //loop over quadrature points in element
        for(unsigned qpi=0; qpi<e->n_quadPoints; ++qpi) {
            FullMatrix shapeGrads = e->grads[qpi];
            Vector shapeVals = e->vals[qpi];
            double jxw =  e->JxW(qpi);

            //assemble element tangent matrix
            for(int A=0; A<Ndof; ++A) {
                for(int B=A; B<Ndof; ++B) {
                    for(unsigned i=0; i<e->n_spaceDim; ++i) {
                        Kel(A,B) += shapeGrads(A,i)*shapeGrads(B,i)*jxw;
                    }
                    Kel(B,A) = Kel(A,B); //exploit symmetry
                }

                //create element force vector
                Fel(A) += fvol*shapeVals(A)*jxw;
            }

        } // end quadrature point loop
    
    
        //assemble element Tangent/force into global
        for(int A=0; A<Ndof; ++A) {
            int GA = e->global_dofs[A];
            for(int B=0; B<Ndof; ++B) {
                int GB = e->global_dofs[B];

                Kcoo.add(GA, GB, Kel(A,B));
            }
      
            Fsys(GA) += Fel(A);
        }
    
    } // end element loop
  
  
    std::cout << "Assembly done\n";

}


void Poisson::solve() {

    sparse_csr Ksys(Kcoo);
    /*
      Vector u_bc(M.get_n_nodes(), 0.0);
      for(unsigned i=0; i<bcnodes.size(); ++i) {
      u_bc(bcnodes[i]) = bcvals[i];
      }

      Vector f_bc = Ksys*u_bc;
      Fsys += f_bc*(-1);
  
      for(unsigned i=0; i<bcnodes.size(); ++i) {
      int dof = bcnodes[i];
      Ksys.zero_row(dof);
      Ksys.zero_col(dof);
      Ksys.assign(dof, dof, 1.0);
      Fsys(dof) = bcvals[i];
      }
    */
    BC.apply(Ksys, Fsys);
    BC2.apply(Ksys,Fsys);
    Vector u_bc = BC.getUbc() + BC2.getUbc();
  
    std::cout << "BCs applied. Solving...\n";
  
    // apply linear solver (conjugate gradient)
    //ILUPreconditioner Kilu(Ksys);
    //Usol = pcg_solve(Ksys, Fsys, u_bc,Kilu);
    Usol = cg_solve(Ksys,Fsys,u_bc);
    Usol.print();
    std::cout << "Solved!\n";
}

void Poisson::output(const std::string & outputFilename) {

    VTKMesh vtkm(EF);
    VTKScalarData vtku(Usol, VTKObject::VTKPOINTDATA, std::string("Usol"));
  
    VTKOutput vo;
    vo.addVTKObject(&vtkm);
    vo.addVTKObject(&vtku);
  
    vo.write(outputFilename);
}

void Poisson::run(const std::string & outputFilename) {
  
    setup();//inputFilename);
    assemble();
    solve();
    output(outputFilename);

}


YAFEL_NAMESPACE_CLOSE
