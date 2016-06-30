#ifndef _YAFEL_GENERIC_MESH_HPP
#define _YAFEL_GENERIC_MESH_HPP

/*
 * GenericMesh: Base class utilizing CRTP to define the interface for finite element
 * mesh data structures. It cannot be instantiated directly, but enables static
 * polymorphism for all functions taking a reference to a mesh.
 */

#include "yafel_globals.hpp"
#include "lin_alg/tensor/Tensor.hpp"
#include "mesh/MeshIterator.hpp"
#include "utils/ElementType.hpp"
#include "mesh/Face.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <exception>

YAFEL_NAMESPACE_OPEN

template<typename T, unsigned NSD>
class GenericMesh {

public:
    using coordinate_type = Tensor<NSD,1,double>;
    //using size_type = typename yafel::size_type;
    using element_container = std::vector<size_type>;

    // vector of all Faces in mesh. Only constructed if build_faces() is called.
    std::vector<Face> mesh_faces;

    //pointers to faces for each cell (offsets in mesh_faces vector)
    std::vector<std::vector<size_type> > cell_faces; 

    //boolean to indicate if build_faces has been called
    bool faces_built;


    // Functions that all meshes must support
    inline size_type n_nodes() const {return static_cast<T const&>(*this).n_nodes();}
    inline size_type n_elements() const {return static_cast<T const&>(*this).n_elements();}
    inline ElementType element_type(size_type elnum) const {
	return static_cast<T const&>(*this).element_type(elnum);
    }

    //constructor
    GenericMesh() : mesh_faces(), cell_faces(), faces_built(false) {}
  
    coordinate_type node(size_type nodenum) const {return static_cast<T const&>(*this).node(nodenum);}
    element_container element(size_type elnum) const {return static_cast<T const&>(*this).element(elnum);}

    /*
     * Function that is user-called if they want the internal faces to be built
     * and stored in internal_faces.
     */
    void build_faces();
  
    //utility function to return number of faces for a given cell
    size_type n_cell_faces(size_type elnum);

    //return nodes blonging to fnum-th face of element elnum (ElementType-dependent)
    element_container face_nodes(size_type elnum, size_type fnum);


    // iterator begin() and end(). This should not be implemented in children of this class.
    MeshIterator<T> begin() const {
	return MeshIterator<T>(*this, 0);
    }
    MeshIterator<T> end() const {
	return MeshIterator<T>(*this, n_elements());
    }
  
  
    // Stuff for CRTP to work
    operator T&(){return static_cast<T&>(*this);}
    operator T const&() const {return static_cast<T const&>(*this);}
};


//=================================================================
/*
 * Implementation
 */
//=================================================================

template<typename MT, unsigned NSD>
void GenericMesh<MT,NSD>::build_faces() {
  
  
    //filter elements for only volume-integrable elements.
    // ie, in 2d, filter for only quads/tris,
    // in 3d, filter for tets, hexes
    std::vector<size_type> vol_elements;
    cell_faces.resize(n_elements());
    for(size_type e=0; e<n_elements(); ++e) {
	ElementType et = element_type(e);
	switch(et) {
	case ElementType::LINEAR_LINE:
	case ElementType::QUADRATIC_LINE:
	case ElementType::CUBIC_LINE:
	    if(NSD==1) {
		vol_elements.push_back(e);
		cell_faces[e].resize(n_cell_faces(e),0);
		break;
	    }
	case ElementType::LINEAR_TRI:
	case ElementType::QUADRATIC_TRI:
	case ElementType::CUBIC_TRI:
	case ElementType::LINEAR_QUAD:
	case ElementType::QUADRATIC_QUAD:
	case ElementType::CUBIC_QUAD:
	    if(NSD==2) {
		vol_elements.push_back(e);
		cell_faces[e].resize(n_cell_faces(e),0);
		break;
	    }
	case ElementType::LINEAR_TET:
	case ElementType::QUADRATIC_TET:
	case ElementType::CUBIC_TET:
	case ElementType::LINEAR_HEX:
	case ElementType::QUADRATIC_HEX:
	case ElementType::CUBIC_HEX:
	    if(NSD==3) {
		vol_elements.push_back(e);
		cell_faces[e].resize(n_cell_faces(e),0);
		break;
	    }
	case ElementType::NULL_ELEMENT:
	default:
	    break;
	}
    }
  
  
    //loop over elements and build node-element adjacency sets
    std::vector<std::set<size_type> > node_sets(n_nodes());
    for(auto e : vol_elements) {
	element_container elem = element(e);

	for(auto n : elem) {
	    node_sets[n].insert(e);
	}
    }


    //loop over elements and build faces
    for(auto e : vol_elements) {
	for(size_type f=0; f<n_cell_faces(e); ++f) {
	    auto f_nodes = face_nodes(e,f);
	    auto nset = node_sets[f_nodes[0]];
	    std::vector<size_type> res(nset.size());
	    std::copy(nset.begin(), nset.end(), res.begin());
	    for(size_type i=1; i<f_nodes.size(); ++i) {
		auto endit = std::set_intersection(res.begin(), res.end(),
						   node_sets[f_nodes[i]].begin(), node_sets[f_nodes[i]].end(),
						   res.begin());
		res.resize(endit-res.begin());
	    }
      
	    if(res.size() == 1) {
		//boundary face
		Face F(f_nodes, e, f, e, f, true);
		size_type fi_global = mesh_faces.size();
		mesh_faces.push_back(F);
		cell_faces[e][f] = fi_global;
	    }
	    else if(res.size() == 2) {
		//construct face
		size_type adj; //, tmp1, tmp2;
		adj = (e==res[0]) ? res[1] : res[0];

		//if e > adj, do nothing
		if(e > adj) {
		    continue;
		}


		//find local face number of F in adjacent cell
		std::sort(f_nodes.begin(), f_nodes.end());
		for(size_type af=0; af<n_cell_faces(adj); ++af) {
		    auto af_nodes = face_nodes(adj,af);
		    std::sort(af_nodes.begin(), af_nodes.end());
		    bool sameface = true;
		    for(size_type i=0; i<f_nodes.size(); ++i) {
			sameface = sameface && f_nodes[i]==af_nodes[i];
		    }

		    if(sameface) {

			Face F(f_nodes, e, f, adj, af, false);
			size_type fi_global = mesh_faces.size();
			mesh_faces.push_back(F);
			cell_faces[e][f] = fi_global;
			cell_faces[adj][af] = fi_global;
			break;
		    }
		}
        
	    }
	    else {
		throw std::runtime_error("More than 2 elements on face. Something is wrong.");
	    }
	}

    }

    //set the "faces_built" flag so this doesn't have to be done again
    faces_built = true;
} //end build_faces()
  
  
  
//==================================================================
template<typename MT, unsigned NSD>
size_type 
GenericMesh<MT,NSD>::n_cell_faces(size_type elnum) {
    ElementType et = element_type(elnum);
    switch(et) {
    case ElementType::LINEAR_LINE:
    case ElementType::QUADRATIC_LINE:
    case ElementType::CUBIC_LINE:
	return 2;
    case ElementType::LINEAR_TRI:
    case ElementType::QUADRATIC_TRI:
    case ElementType::CUBIC_TRI:
	return 3;
    case ElementType::LINEAR_QUAD:
    case ElementType::QUADRATIC_QUAD:
    case ElementType::CUBIC_QUAD:
    case ElementType::LINEAR_TET:
    case ElementType::QUADRATIC_TET:
    case ElementType::CUBIC_TET:
	return 4;
    case ElementType::LINEAR_HEX:
    case ElementType::QUADRATIC_HEX:
    case ElementType::CUBIC_HEX:
	return 6;
    case ElementType::NULL_ELEMENT:
    default:
	return 0;
    }

    return 0;
}

//==================================================================
template<typename MT, unsigned NSD>
typename GenericMesh<MT,NSD>::element_container
GenericMesh<MT,NSD>::face_nodes(size_type elnum, size_type fnum) {

    auto e = element(elnum);
    switch(element_type(elnum)) {
    case ElementType::LINEAR_LINE:
    case ElementType::QUADRATIC_LINE:
    case ElementType::CUBIC_LINE:
	return element_container{fnum};
    case ElementType::LINEAR_TRI:
    case ElementType::QUADRATIC_TRI:
    case ElementType::CUBIC_TRI:
	return element_container{e[fnum], e[(fnum+1)%3]};
    case ElementType::LINEAR_QUAD:
    case ElementType::QUADRATIC_QUAD:
    case ElementType::CUBIC_QUAD:
	return element_container{e[fnum], e[(fnum+1)%4]};
    
	//not handling 3D yet
    case ElementType::LINEAR_TET:
    case ElementType::QUADRATIC_TET:
    case ElementType::CUBIC_TET:
    case ElementType::LINEAR_HEX:
    case ElementType::QUADRATIC_HEX:
    case ElementType::CUBIC_HEX:
    case ElementType::NULL_ELEMENT:
    default:
	throw std::invalid_argument("GenericMesh::face_nodes: 3D Elements currently not supported");
    }

    return element_container{};
}

YAFEL_NAMESPACE_CLOSE

#endif
