#include "yafel.hpp"
#include <iostream>
#include <algorithm>
#include <vector>

using namespace yafel;

constexpr unsigned NSD=3;

template<typename MT>
void color_mesh(const GenericMesh<MT,NSD> &M) {

  std::vector< std::vector<std::size_t> > node_elem_adjacency(M.n_nodes());
  //build adjacency list
  for(std::size_t i=0; i<M.n_elements(); ++i) {

    auto elem = M.element(i);

    for(auto n : elem) {
      node_elem_adjacency[n].push_back(i);
    }
  }

  std::vector<double> colors(M.n_elements(), -1);
  
  //loop over elements and greedily color
  for(std::size_t elnum=0; elnum < M.n_elements(); ++elnum) {

    std::vector<double> used_colors;
    used_colors.reserve(20);
    
    for(auto n : M.element(elnum)) {
      for(auto e : node_elem_adjacency[n]) {
	auto c = colors[e];
	if(c != -1) {
	  used_colors.push_back(c);
	}
      }
    }
    std::sort(begin(used_colors), end(used_colors));
    auto newend = std::unique(begin(used_colors), end(used_colors));

    double c = -1;
    bool color_set = false;
    std::size_t i=0;
    for(auto it=begin(used_colors); it != newend; ++it, ++i) {
      c = static_cast<double>(i);
      if(c < used_colors[i]) {
	colors[elnum] = static_cast<double>(i);
	color_set = true;
	break;
      }
    }
    if(!color_set) {
      colors[elnum] = c+1;
    }
    
  }


  //output
  VTKOutput vout;
  VTKMesh<MT,NSD> vtkm(M);
  vout.addVTKObject(&vtkm);
  VTKScalarData dat(colors, VTKObject::VTKCELLDATA, "color");
  vout.addVTKObject(&dat);

  vout.write("meshcolors.vtu");
}


int main(int argc, char **argv) {

  if(argc < 2) {
    std::cout << "Mesh file not specified, using 10x10 RectilinearMesh<2>.\n";
    RectilinearMesh<NSD> M(std::vector<double>(NSD,1),
			   std::vector<std::size_t>(NSD,10));

    color_mesh(M);

    return 0;
  }
  else {
    // Read a gmsh mesh
    GmshMesh<NSD> M(argv[1]);

    color_mesh(M);
  }

  return 0;
}
