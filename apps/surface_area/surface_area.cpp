#include "yafel.hpp"
#include <iostream>

using namespace yafel;

constexpr int NSD = 3;

int main() {

    GmshMesh<NSD> M("cube.msh");
    
    CG_DoFManager<GmshMesh<NSD>,NSD> dofm(M,1);

    ElementFactory<GmshMesh<NSD>,NSD> EF(M,dofm);

    double SurfaceArea = 0;
    for(std::size_t e=0; e<M.n_elements(); ++e) {

	auto &E = EF.getElement(e);
	if(E.n_topoDim != NSD-1) {
	    continue;
	}
	std::cout << e << std::endl;
	
	E.update_element(M,e);


	for(std::size_t qpi=0; qpi<E.n_quadPoints; ++qpi) {
	    SurfaceArea += 1.0*E.JxW(qpi);
	}
    }

    std::cout << "Area = " << SurfaceArea << std::endl;

    return 0;
}
