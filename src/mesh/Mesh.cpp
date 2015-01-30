#include "mesh/Mesh.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <queue>
#include <utility>

YAFEL_NAMESPACE_OPEN
//==================================================================================
Mesh::Mesh() : n_elems(0), n_nodes(0), minLength(0) {}

//==================================================================================
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

//==================================================================================
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

//==================================================================================
// This next function is my first stab at using C++ iterators effectively. I use them
// where possible, but occasionally I need the numeric value of an array index, so
// "unsigned" vars must be used.
//==================================================================================
void Mesh::reorder_rcm() {
  
  // make node adjacency list
  std::vector<std::vector<unsigned> > adj_list(get_n_nodes());
  std::vector<unsigned> degree(get_n_nodes(), 0);
  std::vector<bool> visited(get_n_nodes(), false);
  
  std::vector<unsigned> R; //results vector R[i] is the old node number for new node i
  std::queue<unsigned> Q; // fifo queue for breadth-first search
  
  for(auto el = elements.begin(); el<elements.end(); ++el) {
    for(unsigned i=0; i<(*el).size(); ++i) {
      for(unsigned j=i+1; j<(*el).size(); ++j) {
	adj_list[(*el)[i]].push_back((*el)[j]);
	adj_list[(*el)[j]].push_back((*el)[i]);
      }
    }
  }
  
  //compute degrees of each node, locate max-degree node (imax)
  auto deg_it = degree.begin();
  unsigned imax = 0;
  unsigned max_deg = 0;
  unsigned i=0;
  for(auto it=adj_list.begin(); it<adj_list.end(); ++it, ++deg_it, ++i) {
    std::sort((*it).begin(), (*it).end());
    auto end = std::unique((*it).begin(), (*it).end());
    *deg_it = end-(*it).begin();
    max_deg = (max_deg > *deg_it) ? max_deg : *deg_it;
    imax = (max_deg > *deg_it) ? imax : i;
  }

  //breadth-first search, adding high-degree neighbors first
  Q.push(imax);
  
  while(!Q.empty()) {
    unsigned node = Q.front();
    Q.pop();
    if(visited[node])
      continue;

    R.push_back(node);
    visited[node] = true;
    
    std::vector<std::pair<unsigned, unsigned> > neighbors;
    for(unsigned i=0; i<degree[node]; ++i) {
      if(i != node) {
	neighbors.push_back(std::pair<unsigned,unsigned>(degree[adj_list[node][i]], adj_list[node][i]));
      }
    }
    std::sort(neighbors.begin(), neighbors.end());

    //append indices to queue in descending order of degree (std::sort puts them in ascending)
    for(auto it=neighbors.rbegin(); it<neighbors.rend(); ++it) {
      Q.push(it->second);
    }
  }
  
  // create new nodes vector and update elements 
  std::vector<Vector> new_nodes;
  for(unsigned n=0; n<get_n_nodes(); ++n) {
    new_nodes.push_back(nodal_coords[R[n]]);
  }
  nodal_coords = new_nodes;

  std::vector<unsigned> Rinv(get_n_nodes(),0); //maps: Rinv[oldNodeNum] -> newNodeNum
  for(unsigned i=0; i<get_n_nodes(); ++i) {
    Rinv[R[i]] = i;
  }

  for(auto e=elements.begin(); e<elements.end(); ++e) {
    for(auto n=e->begin(); n<e->end(); ++n) {
      *n = Rinv[*n];
    }
  }
  
}


//==================================================================================
//==================================================================================
//==================================================================================
YAFEL_NAMESPACE_CLOSE
