#include "mesh/Mesh.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>

YAFEL_NAMESPACE_OPEN

Mesh::Mesh() : n_elems(0), n_nodes(0), minLength(0) {}

Mesh::Mesh(const std::vector<Vector> & nodes,
	   const std::vector< std::vector<int> > & elems,
	   const std::vector<int> & eltype) {
  this->nodal_coords = nodes;
  this->elements = elems;
  this->element_type = eltype;

  this->n_elems = elems.size();
  this->n_nodes = nodes.size();
}

Mesh::Mesh(const std::vector< Vector > & nodes,
	   const std::vector< std::vector<int> > & elems,
	   const std::vector<int> & eltype,
	   const std::vector< std::vector<int> > & _tags) {
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
  minLength = std::sqrt(nodal_coords[0].dot(nodal_coords[1]));
  for(unsigned elem=0; elem<n_elems; ++elem) {
    
    for(unsigned i=0; i<elements[elem].size(); ++i) {
      for(unsigned j=0; j<elements[elem].size(); ++j) {

	double d = std::sqrt(nodal_coords[elements[elem][i]].dot(nodal_coords[elements[elem][j]]));
	
	minLength = (minLength<=d) ? minLength : d;
      }
    }
    
  }//end element loop
  
}



YAFEL_NAMESPACE_CLOSE
