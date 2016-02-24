#ifndef _YAFEL_VTKVECTORDATA_HPP
#define _YAFEL_VTKVECTORDATA_HPP

/*
 * Data type representing the output of a vector field. It is intended
 * to be used for vectors of the same dimension as the spatial dimension.
 *
 * Due to constraints with the VTK output, it can only be used for up to NSD=3
 */

#include "yafel_globals.hpp"
#include "output/VTKObject.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class VTKVectorData : public VTKObject {

private:
  std::vector<Tensor<NSD,1> > data;

public:
  VTKVectorData(const std::vector<Tensor<NSD,1> > &d, VTKObject::VTKObjectType ot, const std::string &s);

  void write(FILE *fp);
};



/*
 * Implementation
 */ 
template<unsigned NSD>
VTKVectorData<NSD>::VTKVectorData(const std::vector<Tensor<NSD,1> > &d, 
                                    VTKObject::VTKObjectType ot, const std::string &s) :
                                  VTKObject(VTKObject::VTKVECTORDATA, ot, s), data(d) {}


template<unsigned NSD>
void VTKVectorData<NSD>::write(FILE *fp) {
  
  fprintf(fp, "<DataArray type=\"Float32\" format=\"ascii\" Name=\"%s\" NumberOfComponents=\"3\">\n",
	  getName().c_str());
  
  for(unsigned i=0; i<data.size(); ++i) {
    
    fprintf(fp, "%f %f %f\n", 
	    data[i](0),
	    (NSD>=2)?data[i](1):0,
	    (NSD>=3)?data[i](2):0);
  }

  fprintf(fp, "</DataArray>\n");
}


YAFEL_NAMESPACE_CLOSE

#endif
