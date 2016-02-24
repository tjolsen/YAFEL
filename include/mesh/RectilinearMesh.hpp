#ifndef __YAFEL_RECTILINEAR_MESH
#define __YAFEL_RECTILINEAR_MESH

/*
 * RectilinearMesh:
 *
 * Data structure for N-dimensional rectilinear meshes.
 * Users are able to specify the dimensions and mesh resolution
 * along each dimension.
 * The mesh is assumed to be aligned with coordinate axes, with
 * nodes[0] corresponding to X = (0,0,...).
 *
 * Nodes and elements are implicitly stored. Values for a node
 * coordinate or element connectivity are created on demand and
 * returned by value.
 *
 *
 * Elements are reported using the traditional connectivity,
 * (ie, counterclockwise for 2D, two ccw layers for 3D).
 * The same connectivity is used for equivalent gmsh meshes.
 *
 *
 * Currently, elements are assumed to be linear, but this should
 * be updated.
 */

#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/ElementType.hpp"

#include <vector>
#include <exception>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class RectilinearMesh : public GenericMesh<RectilinearMesh<NSD>, NSD> {
  
public:
  using coordinate_type = typename GenericMesh<RectilinearMesh<NSD>,NSD>::coordinate_type;
  using size_type = typename GenericMesh<RectilinearMesh<NSD>,NSD>::size_type;
  using element_container = typename GenericMesh<RectilinearMesh<NSD>,NSD>::element_container;
  
  // mesh characteristics
  std::vector<double> dims; //hyper-rectangle dimensions
  std::vector<size_type> dim_elems; // <-- number of elements along each dimension; dx = dims[i]/dim_elems[i]
  

  // default is 1 x 1 x ... cube with 10 elements along each dimension
  RectilinearMesh(const std::vector<double> & _dims,
                  const std::vector<size_type> & _dim_elems) :
    dims(_dims), dim_elems(_dim_elems)
  {
    // only supporting 1D, 2D, 3D for now. Unlikely to need anything else.
    static_assert(NSD>0 && NSD<4, "Rectilinear mesh currently intended for NSD={1,2,3}");
    
    if(dims.size() != NSD) {
      throw(std::invalid_argument("RectilinearMesh: wrong size _dims argument"));
    }
    if(dim_elems.size() != NSD) {
      throw(std::invalid_argument("RectilinearMesh: wrong size _dim_elems argument"));
    }
  }

  RectilinearMesh() : RectilinearMesh(std::vector<double>(NSD,1.0),
                                      std::vector<size_type>(NSD,10))
  {}


  // Interface implementation
  
  inline size_type n_nodes() const {
    size_type ret(1);
    for(auto de : dim_elems) {
      ret *= (de+1);
    }
    return ret;
  }

  inline size_type n_elements() const {
    size_type ret(1);
    for(auto de : dim_elems) {
      ret *= de;
    }
    return ret;
  }

  inline ElementType element_type(size_type) const {return ElementType::LINEAR_QUAD;}
  
  inline coordinate_type node(size_type nodenum) const {
    
    coordinate_type ret;
    size_type stride(1);
    for(size_type i=0; i<NSD; ++i) {
      //get mesh spacing
      double dx = dims[i]/dim_elems[i];
      
      // get node index along dimension i
      size_type idx = (nodenum/stride) % (dim_elems[i]+1);
      
      //set coordinate
      ret(i) = dx*idx;

      // update stride for next dim
      stride *= (dim_elems[i]+1);
    }
    
    return ret;
  }
  
  inline element_container element(size_type elnum) const;
  
private:
  
  
};


/*
 * Implementing elements() function separately for 1D, 2D, and 3D, since it would
 * be complicated to form the correct element connectivity in a general way
 */
template<>
inline typename RectilinearMesh<1>::element_container
RectilinearMesh<1>::element(size_type elnum) const {
  return element_container{elnum, elnum+1};
}

template<>
inline typename RectilinearMesh<2>::element_container
RectilinearMesh<2>::element(size_type elnum) const {
  
  size_type el_x = elnum % dim_elems[0];
  size_type el_y = (elnum/dim_elems[0]) % dim_elems[1];
  
  //id of lower-left node
  size_type LL = el_y*(dim_elems[0]+1) + el_x;

  return element_container{LL, LL+1, LL+dim_elems[0]+2, LL+dim_elems[0]+1};
}


template<>
inline typename RectilinearMesh<3>::element_container
RectilinearMesh<3>::element(size_type elnum) const {

  size_type y_stride = dim_elems[0]+1;
  size_type z_stride = (dim_elems[0]+1)*(dim_elems[1]+1);
  
  size_type el_x = elnum % dim_elems[0];
  size_type el_y = (elnum/dim_elems[0]) % dim_elems[1];
  size_type el_z = (elnum/dim_elems[0]/dim_elems[1]) % dim_elems[2];

  //id of lower-left node
  size_type LL = el_z*z_stride + el_y*y_stride + el_x;

  return element_container{LL, LL+1, 
      LL+y_stride+1, LL+y_stride,
      LL+z_stride, LL+z_stride+1, 
      LL+z_stride+y_stride+1, LL+z_stride+y_stride};
}



YAFEL_NAMESPACE_CLOSE

#endif
