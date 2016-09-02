#ifndef _YAFEL_GMSHMESH_HPP
#define _YAFEL_GMSHMESH_HPP


/* 
 * GmshMesh:
 * 
 * This data structure represents a mesh read from a gmsh .msh file.
 * Since these meshes have heterogeneous element types, it must explicitly
 * store the node locations, element connectivity, element types, and "tags"
 * (see gmsh documentation).
 * 
 */

#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include "utils/ElementType.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <string.h>

YAFEL_NAMESPACE_OPEN

template<unsigned NSD>
class GmshMesh : public GenericMesh<GmshMesh<NSD>, NSD> {

public:
    using coordinate_type = typename GenericMesh<GmshMesh<NSD>,NSD>::coordinate_type;
    using element_container = typename GenericMesh<GmshMesh<NSD>,NSD>::element_container;

    std::vector<coordinate_type> _nodes;
    std::vector<element_container> _elements;
    std::vector<size_type> _element_type;
    std::vector<std::vector<size_type> > _element_tags;

    // Constructors
    GmshMesh(const std::string &fname) {
    
	std::ifstream in(fname.c_str());
	if(!in.good())
	    throw(std::runtime_error("GmshMesh: Could not open mesh file"));

	load_from_stream(in);

	in.close();
    }


    GmshMesh(std::istream &in) {
	load_from_stream(in);
    }

    // Interface implementation
    inline size_type n_nodes() const {
	return _nodes.size();
    }
  
    inline size_type n_elements() const {
	return _elements.size();
    }

    inline coordinate_type node(size_type nodenum) const {
	return _nodes[nodenum];
    }

    inline element_container element(size_type elnum) const {
	return _elements[elnum];
    }

    inline ElementType element_type(size_type elnum) const {
	return ElementType_Mappings::gmsh_to_ElementType(_element_type[elnum]);
    }

private:
    void load_from_stream(std::istream &in);

    // little utility function to split strings and store results into a vector
    void string_split(const std::string &s, char delim, 
		      std::vector<std::string> & elems) {
	std::stringstream ss(s);
	std::string word;
	while(std::getline(ss, word, delim)) {
	    elems.push_back(word);
	}
    
	return;
    }
  
  
}; //end class



/*
 * load_from_stream() implementation
 */
template<unsigned NSD>
void GmshMesh<NSD>::load_from_stream(std::istream &in) {
  
    if(!in.good()) {
	throw(std::runtime_error("GmshMesh: bad std::istream"));
    }
  
    typedef enum action {
	ACTION_UNSET,
	PARSE_NODES,
	PARSE_ELEMENTS
    } action_t;
  
  
    //  std::vector<coordinate_type> nodal_coords;
    //  std::vector<std::vector<size_type> > elements;
    //  std::vector<size_type> el_type;
    //  std::vector<std::vector<size_type> > tags;
  
    action_t currentAction = ACTION_UNSET;

    while(!in.eof()) {
	std::string line;
	std::getline(in, line);
	if(in.eof()) {
	    break;
	}
	std::vector<std::string> words;
	string_split(line, ' ', words);
    
	if( strcmp(words.at(0).c_str(), "$Nodes")==0 ) {
	    currentAction = PARSE_NODES;
	    std::getline(in, line);
	    int nNodes = atoi(line.c_str());
	    _nodes.resize(nNodes);
	    continue;
	} 
	else if( strcmp(words.at(0).c_str(), "$EndNodes")==0 ) {
	    currentAction = ACTION_UNSET;
	    continue;
	}
	else if( strcmp(words.at(0).c_str(), "$Elements")==0 ) {
	    currentAction = PARSE_ELEMENTS;
	    std::getline(in, line);
	    int nElems = atoi(line.c_str());
	    _elements.resize(nElems);
	    _element_tags.resize(nElems);
	    _element_type.resize(nElems);
	    continue;      
	}
	else if( strcmp(words.at(0).c_str(), "$EndElements")==0 ) {
	    currentAction = ACTION_UNSET;
	    continue;
	}
    
	int id = -1;
	switch (currentAction) {
	case ACTION_UNSET:
	    break;
      
	case PARSE_NODES:
	{
	    id = atoi(words[0].c_str());
	    coordinate_type node;
	    for(unsigned i=0; i<NSD; ++i) {
		node(i) = atof(words[i+1].c_str());
	    }
	    _nodes[id-1] = node;
	    break;
	}
	case PARSE_ELEMENTS:
	{
	    id = atoi(words[0].c_str());
	    _element_type[id-1] = atoi(words[1].c_str());
	    int ntags = atoi(words[2].c_str());
	    int nodes_in_el = words.size() - ntags - 3;
	    std::vector<size_type> tag(ntags,0);
	    for(int i=3; i<3+ntags; ++i) {
		tag[i-3] = atoi(words[i].c_str());
	    }
	    _element_tags[id-1] = tag;
	    element_container el(nodes_in_el, 0); //init to 0. gonna fail hard somewhere, if it does.
	    for(size_type i=3+ntags; i<words.size(); ++i) {
		el[i-ntags-3] = atoi(words[i].c_str())-1; //using 0-based node numbering from gmsh 1-based
	    }
	    _elements[id-1] = el;
	    break;
	}
	}
    
    } //end while(!in.eof())
  
  
}



YAFEL_NAMESPACE_CLOSE

#endif
