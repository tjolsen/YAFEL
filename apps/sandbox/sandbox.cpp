#include "yafel.hpp"
#include <vector>
#include <iostream>

using namespace yafel;

constexpr int NSD = 3;

Tensor<3,1> cross(Tensor<3,1> lhs, Tensor<3,1> rhs) {

    return Tensor<3,1>{ lhs(1)*rhs(2)-lhs(2)*rhs(1),
	    -(lhs(0)*rhs(2) - rhs(0)*lhs(2)),
	    lhs(0)*rhs(1) - rhs(0)*lhs(1)
	    };
    
}


void func() {

    
    GmshMesh<NSD> M("square.msh");
    

    CG_DoFManager<GmshMesh<NSD>,NSD> dofm(M,1);

    ElementFactory<GmshMesh<NSD>,NSD> EF(M,dofm);

    for(std::size_t i=0; i<M.n_elements(); ++i) {
	auto &E = EF.getElement(i);
	if(E.n_topoDim != 3) {
	    continue;
	}
	E.update_element(M,i);

	for(auto n : M.element(i)) {
	    for(auto x : M.node(n)) {
		std::cout << x << "  ";
	    }
	    std::cout << std::endl;
	}

	
	Tensor<NSD,1> grad;
	//Tensor<NSD,1>xi{1.0/2.0,1.0/2.0,1.0/2.0};
 	Tensor<NSD,1> xi{ 1.0/3.0 , 0, 1.0/3.0  };
	
	//xi at center of angled
	auto J = E.calcJ_xi(xi);
	for(auto j : J)
	    std::cout << j << std::endl;
	
	Tensor<NSD,1> a{J(0,1), J(1,1), J(2,1)};
	Tensor<NSD,1> b{J(0,2), J(1,2), J(2,2)};
	
	std::cout << std::endl;
	for(auto v : cross(a,b)) {
	    std::cout << v << "  ";
	}
	std::cout << std::endl << "-----------------------------------------" << std::endl;
	break;
    }
}



int main() {

    RectilinearMesh<NSD> M(std::vector<double>(NSD,1),
			   std::vector<std::size_t>(NSD,1));


    CG_DoFManager<RectilinearMesh<NSD>,NSD> dofm(M,1);

    ElementFactory<RectilinearMesh<NSD>,NSD> EF(M,dofm);

    auto & E = EF.getElement(0);

    E.update_element(M,0);

    Tensor<NSD,1> grad,xi{1,0,0};

    //xi at center of right face
    auto J = E.calcJ_xi(xi);
    for(auto j : J)
	std::cout << j << std::endl;

    Tensor<NSD,1> a{J(0,1), J(1,1), J(2,1)};
    Tensor<NSD,1> b{J(0,2), J(1,2), J(2,2)};

    std::cout << std::endl;
    for(auto v : cross(a,b)) {
	std::cout << v << "  ";
    }
    std::cout << std::endl << "-----------------------------------------" << std::endl;

    


    func();

    return 0;
}
