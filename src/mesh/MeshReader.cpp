#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "mesh/MeshReader.hpp"
#include "mesh/Mesh.hpp"

YAFEL_NAMESPACE_OPEN

void MeshReader::string_split(std::string s, char delim, 
			     std::vector<std::string> & elems) {
  std::stringstream ss(s);
  std::string word;
  while(std::getline(ss, word, delim)) {
    elems.push_back(word);
  }
  
  return;
}

Mesh MeshReader::gmsh_read(std::string fname) {
  
  std::ifstream infile(fname.c_str());
  if(!infile.is_open()) {
    std::cerr << "Could not open gmsh .msh file" << std::endl;
  }
  
  typedef enum action {
    ACTION_UNSET,
    PARSE_NODES,
    PARSE_ELEMENTS
  } action_t;
  
  
  std::vector< Vector > nodal_coords;
  std::vector<std::vector<unsigned> > elements;
  std::vector<unsigned> el_type;
  std::vector<std::vector<unsigned> > tags;

  action_t currentAction = ACTION_UNSET;
  
  while(!infile.eof()) {
    std::string line;
    std::getline(infile, line);
    if(infile.eof()) {
      break;
    }
    std::vector<std::string> words;
    string_split(line, ' ', words);
    
    if( strcmp(words.at(0).c_str(), "$Nodes")==0 ) {
      currentAction = PARSE_NODES;
      std::getline(infile, line);
      int nNodes = atoi(line.c_str());
      nodal_coords.resize(nNodes);
      continue;
    } 
    else if( strcmp(words.at(0).c_str(), "$EndNodes")==0 ) {
      currentAction = ACTION_UNSET;
      continue;
    }
    else if( strcmp(words.at(0).c_str(), "$Elements")==0 ) {
      currentAction = PARSE_ELEMENTS;
      std::getline(infile, line);
      int nElems = atoi(line.c_str());
      elements.resize(nElems);
      tags.resize(nElems);
      el_type.resize(nElems);
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
	Vector node(3,0.0);
	for(int i=0; i<3; ++i) {
	  node(i) = atof(words[i+1].c_str());
	}
	nodal_coords[id-1] = node;
	break;
      }
    case PARSE_ELEMENTS:
      {
	id = atoi(words[0].c_str());
	el_type[id-1] = atoi(words[1].c_str());
	int ntags = atoi(words[2].c_str());
	int nodes_in_el = words.size() - ntags - 3;
	std::vector<unsigned> tag(ntags,0);
	for(int i=3; i<3+ntags; ++i) {
	  tag[i-3] = atoi(words[i].c_str());
	}
	tags[id-1] = tag;
	std::vector<unsigned> el(nodes_in_el, -1); //init to -1 so fails hard,if at all
	for(unsigned i=3+ntags; i<words.size(); ++i) {
	  el[i-ntags-3] = atoi(words[i].c_str())-1; //using 0-based node numbering
	}
	elements[id-1] = el;
	break;
      }
    }
    
  } //end while(!infile.eof())

  Mesh M(nodal_coords,elements, el_type, tags);
  return M;
}



YAFEL_NAMESPACE_CLOSE
