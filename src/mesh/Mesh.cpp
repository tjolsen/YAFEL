#include "mesh/Mesh.hpp"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <queue>
#include <utility>

YAFEL_NAMESPACE_OPEN


//==================================================================================
Mesh::Mesh(const std::vector<coordinate_type> & nodes,
	   const std::vector<std::vector<size_type> > & elems,
	   const std::vector<size_type> & eltype,
	   const std::vector<std::vector<size_type> > & _tags) :
  nodal_coords(nodes), 
  elements(elems), 
  el_tags(_tags), 
  element_type(eltype),
  n_elems(elems.size()), 
  n_nodes(nodes.size())
{}

//==================================================================================
Mesh::Mesh(const std::vector<coordinate_type> & nodes,
	   const std::vector< std::vector<size_type> > & elems,
	   const std::vector<size_type> & eltype) : 
  Mesh(nodes, elems, eltype, std::vector<std::vector<size_type> >(0))
{}

//==================================================================================
Mesh::Mesh() : n_elems(0), n_nodes(0), minLength(0) {}

//==================================================================================
void Mesh::compute_min_length() {

  if(!nodal_coords.size()>=2) {
    perror("Mesh::compute_min_length(): Invalid Mesh. Not enough points");
    exit(1);
  }
  coordinate_type dx0 = (nodal_coords[0]-nodal_coords[1]);
  
  minLength = std::sqrt(contract<1>(dx0, dx0));
  for(size_type elem=0; elem<n_elems; ++elem) {
    if(elements[elem].size() < 2) {
      continue; //point "element"
    }
    
    for(size_type i=0; i<elements[elem].size(); ++i) {
      for(size_type j=i+1; j<elements[elem].size(); ++j) {
	
	coordinate_type dx = nodal_coords[elements[elem][i]] - nodal_coords[elements[elem][j]];
	double d = std::sqrt(contract<1>(dx,dx));
	
	minLength = (minLength <= d) ? minLength : d;
      }
    }
    
  }
  
}

//==================================================================================
// This next function is my first stab at using C++ iterators effectively. I use them
// where possible, but occasionally I need the numeric value of an array index, so
// "size_type" vars must be used.
//==================================================================================
void Mesh::reorder_rcm() {
  
  // make node adjacency list
  std::vector<std::vector<size_type> > adj_list(get_n_nodes());
  std::vector<size_type> degree(get_n_nodes(), 0);
  std::vector<bool> visited(get_n_nodes(), false);
  
  std::vector<size_type> R; //results vector R[i] is the old node number for new node i
  std::queue<size_type> Q; // fifo queue for breadth-first search
  
  for(auto el = elements.begin(); el<elements.end(); ++el) {
    for(size_type i=0; i<(*el).size(); ++i) {
      for(size_type j=i+1; j<(*el).size(); ++j) {
	adj_list[(*el)[i]].push_back((*el)[j]);
	adj_list[(*el)[j]].push_back((*el)[i]);
      }
    }
  }
  
  //compute degrees of each node, locate min-degree node (imin)
  auto deg_it = degree.begin();
  size_type imin = 0;
  size_type min_deg = 1000000; //something unreasonably large so it is overwritten immediately
  size_type i=0;
  for(auto it=adj_list.begin(); it<adj_list.end(); ++it, ++deg_it, ++i) {
    std::sort((*it).begin(), (*it).end());
    auto end = std::unique((*it).begin(), (*it).end());
    *deg_it = end-(*it).begin();
    min_deg = (min_deg < *deg_it) ? min_deg : *deg_it;
    imin = (min_deg < *deg_it) ? imin : i;
  }

  //breadth-first search, adding high-degree neighbors first
  Q.push(imin);
  
  while(!Q.empty()) {
    size_type node = Q.front();
    Q.pop();
    if(visited[node])
      continue;

    R.push_back(node);
    visited[node] = true;
    
    std::vector<std::pair<size_type, size_type> > neighbors;
    for(size_type i=0; i<degree[node]; ++i) {
      if(i != node) {
	neighbors.push_back(std::pair<size_type,size_type>(degree[adj_list[node][i]], adj_list[node][i]));
      }
    }
    std::sort(neighbors.begin(), neighbors.end());

    //append indices to queue in descending order of degree (std::sort puts them in ascending)
    for(auto it=neighbors.begin(); it<neighbors.end(); ++it) {
      Q.push(it->second);
    }
  }
  
  //reverse R
  std::vector<size_type> Rrev(get_n_nodes());
  auto rrev_it = Rrev.begin();
  auto r_it = R.rbegin();
  for(rrev_it = Rrev.begin(); rrev_it < Rrev.end() && r_it<R.rend(); ++rrev_it, ++r_it) {
    *rrev_it = *r_it;
  }
  R.swap(Rrev);

  std::vector<size_type> Rinv(get_n_nodes(),0); //maps: Rinv[oldNodeNum] -> newNodeNum
  for(size_type i=0; i<get_n_nodes(); ++i) {
    Rinv[R[i]] = i;
  }
  
  // create new nodes vector
  std::vector<coordinate_type> new_nodes(get_n_nodes());
  for(size_type n=0; n<get_n_nodes(); ++n) {
    new_nodes[Rinv[n]] = nodal_coords[n];
  }
  nodal_coords.swap(new_nodes);


  //update elements 
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
