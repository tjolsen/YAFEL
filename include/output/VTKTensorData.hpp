#ifndef _YAFEL_VTKTENSORDATA_HPP
#define _YAFEL_VTKTENSORDATA_HPP

#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

// WARNING, VTK CURRENTLY ONLY SUPPORTS 3x3 SYMMETRIC TENSORS

template<unsigned NSD>
class VTKTensorData : public VTKObject {

private:
  std::vector<Tensor<NSD,2> > data;

public:
  VTKTensorData(const std::vector<Tensor<NSD,2> > &d, VTKObject::VTKObjectType ot, const std::string &s);
  
  void write(FILE *fp);
};

/*
 * Implementation
 */
template<unsigned NSD>
VTKTensorData<NSD>::VTKTensorData(const std::vector<Tensor<NSD,2> > &d, 
                                  VTKObject::VTKObjectType ot, const std::string &s) :
  VTKObject(VTKObject::VTKVECTORDATA, ot, s), data(d) {}


template<unsigned NSD>
void VTKTensorData<NSD>::write(FILE *fp) {
  
  fprintf(fp, "<DataArray type=\"Float32\" format=\"ascii\" Name=\"%s\" NumberOfComponents=\"9\">\n",
	  getName().c_str());
  
  for(unsigned i=0; i<data.size(); ++i) {
    
    fprintf(fp, "%f %f %f\n%f %f %f\n%f %f %f\n\n", 
	    (NSD>=1)?data[i](0,0):0, (NSD>=2)?data[i](0,1):0, (NSD>=3)?data[i](0,2):0,
	    (NSD>=2)?data[i](1,0):0, (NSD>=2)?data[i](1,1):0, (NSD>=3)?data[i](1,2):0,
	    (NSD>=3)?data[i](2,0):0, (NSD>=3)?data[i](2,1):0, (NSD>=3)?data[i](2,2):0);
  }
  
  fprintf(fp, "</DataArray>\n");
}


YAFEL_NAMESPACE_CLOSE

#endif
