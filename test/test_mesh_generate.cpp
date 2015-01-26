#include "yafel.hpp"
#include <vector>
#include <string>

using namespace yafel;

int main() {

  MeshGenerator MG(2, std::vector<double>{3,2}, std::vector<unsigned>{51,31});

  Mesh M = MG.getMesh();
  ElementFactory EF(M,DoFManager(1));
  
  VTKOutput vo;
  VTKMesh vtkm(EF);
  vo.addVTKObject(&vtkm);
  
  vo.write(std::string("generated_mesh.vtu"));
  
  return 0;
}
