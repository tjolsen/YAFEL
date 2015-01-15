#include "yafel.hpp"
#include <iostream>
#include <cstring>

using namespace yafel;

int main(int argc, char **argv) {
  
  if(argc < 3) {
    std::cout << "Provide Mesh filename and output basename\n";
    return 1;
  }
  std::string outbase(argv[2]);

  Mesh M(MeshReader::gmsh_read(std::string(argv[1])));
  ElementFactory EF(M,1);

  VTKTimeOutput TO;
  VTKOutput VO;
  VTKMesh vtkm(EF);
  VO.addVTKObject(&vtkm);

  char buf[128];
  for(int t=0; t<10; ++t) {
    sprintf(buf, "%s_%f.vtu", outbase.c_str(), (double)t);
    std::string fname(buf);    

    Vector u(M.get_n_nodes(), t);

    VTKScalarData vtku(u, VTKObject::VTKPOINTDATA, std::string("u"));
    VO.addVTKObject(&vtku);
    VO.write(fname);
    VO.clearData();
    TO.addDataFile(fname, (double)t);    
  }
  std::string time_outbase(outbase);
  time_outbase += ".pvd";
  
  TO.write(time_outbase);
  
  return 0;
}
