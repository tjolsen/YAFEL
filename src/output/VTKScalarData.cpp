#include "output/VTKScalarData.hpp"

YAFEL_NAMESPACE_OPEN

//ctor's
VTKScalarData::VTKScalarData(const std::vector<double> &d, VTKObject::VTKObjectType ot, const std::string &s) :
  VTKObject(VTKObject::VTKSCALARDATA, ot, s), data(d) {}

VTKScalarData::VTKScalarData(const Vector &d, VTKObject::VTKObjectType ot, const std::string &s) :
  VTKObject(VTKObject::VTKSCALARDATA, ot, s)
{
  
  data.resize(d.getLength());
  for(int i=0; i<d.getLength(); ++i) {
    data[i] = d(i);
  }
  
}

void VTKScalarData::write(FILE *fp) {

  fprintf(fp, "<DataArray type=\"Float32\" format=\"ascii\" Name=\"%s\" NumberOfComponents=\"1\">\n",
	  getName().c_str());
  
  for(unsigned i=0; i<data.size(); ++i) {
    fprintf(fp, "%f\n", data[i]);
  }
  
  fprintf(fp, "</DataArray>\n");
  
}

YAFEL_NAMESPACE_CLOSE
