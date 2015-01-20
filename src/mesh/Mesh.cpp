#include "mesh/Mesh.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>

YAFEL_NAMESPACE_OPEN

Mesh::Mesh() : n_elems(0), n_nodes(0), minLength(0) {}

Mesh::Mesh(const std::vector<Vector> & nodes,
	   const std::vector< std::vector<unsigned> > & elems,
	   const std::vector<unsigned> & eltype) {
  this->nodal_coords = nodes;
  this->elements = elems;
  this->element_type = eltype;

  this->n_elems = elems.size();
  this->n_nodes = nodes.size();
}

Mesh::Mesh(const std::vector< Vector > & nodes,
	   const std::vector< std::vector<unsigned> > & elems,
	   const std::vector<unsigned> & eltype,
	   const std::vector< std::vector<unsigned> > & _tags) {
  this->nodal_coords = nodes;
  this->elements = elems;
  this->element_type = eltype;
  this->el_tags = _tags;

  this->n_elems = elems.size();
  this->n_nodes = nodes.size();
}

void Mesh::compute_min_length() {

  if(!nodal_coords.size()>=2) {
    perror("Mesh::compute_min_length(): Invalid Mesh. Not enough points");
    exit(1);
  }
  minLength = (nodal_coords[0]-nodal_coords[1]).norm();
  
  for(unsigned elem=0; elem<n_elems; ++elem) {
    if(elements[elem].size() < 2) {
      continue; //point "element"
    }
    
    for(unsigned i=0; i<elements[elem].size(); ++i) {
      for(unsigned j=i+1; j<elements[elem].size(); ++j) {
	
	Vector dx = nodal_coords[elements[elem][i]] - nodal_coords[elements[elem][j]];
	double d = dx.norm();
	
	minLength = (minLength <= d) ? minLength : d;
      }
    }
    
  }
  
}



YAFEL_NAMESPACE_CLOSE
