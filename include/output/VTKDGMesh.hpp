#ifndef _YAFEL_VTKDGMESH_HPP
#define _YAFEL_VTKDGMESH_HPP

/*
 * VTKDGMesh:
 *
 * Child class of VTKMesh that overrides some methods to allow for
 * a DG-interpretation of a mesh, rather than a CG one.
 * The VTKMesh was modified so that only a few methods need be implemented here.
 *
 * It operates under a fairly restrictive set of assumptions currently:
 *
 *    - It will only work in 2D
 *    - Mesh must consist of quadrilateral elements only. (no lines on boundaries/interfaces)
 *    - Assumes same polynomial order interpolation in all elements
 */

#include "yafel_globals.hpp"
#include "output/VTKMesh.hpp"
#include "utils/ElementVtkType.hpp"
#include "utils/ElementType.hpp"
#include "mesh/GenericMesh.hpp"
#include "element/DG_Quad.hpp"

YAFEL_NAMESPACE_OPEN

template<typename MT, unsigned NSD, typename dataType>
class VTKDGMesh : public VTKMesh<MT,NSD> {

public:
  using size_type = typename VTKMesh<MT,NSD>::size_type;
  using coordinate_type = typename VTKMesh<MT,NSD>::coordinate_type;

  VTKDGMesh(DG_Quad<NSD,dataType> &dgq, GenericMesh<MT,NSD> &M) : VTKMesh<MT,NSD>(M), DGQ(dgq) {}

  //Override VTKMesh virtual functions
  virtual size_type n_nodes() const {return this->M.n_elements()*DGQ.nodes_per_element;}

  virtual coordinate_type node(size_type nodenum) {
    size_type elnum = nodenum/DGQ.nodes_per_element;
    size_type local_nodenum = nodenum % DGQ.nodes_per_element;
    
    DGQ.update_nodes(this->M, elnum);
    
    return DGQ.nodes_x[local_nodenum];
  }
  
  virtual size_type element_size(size_type) const {return DGQ.nodes_per_element;}
  
  virtual size_type element(size_type elnum, size_type locnode) const {
    return elnum*DGQ.nodes_per_element + locnode;
  }

  virtual ElementType element_type(size_type) const {return ElementType::DG_QUAD;}
private:
  DG_Quad<NSD,dataType> & DGQ;
};


YAFEL_NAMESPACE_CLOSE

#endif
