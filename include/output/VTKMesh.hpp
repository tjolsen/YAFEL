#ifndef _YAFEL_VTKMESH_HPP
#define _YAFEL_VTKMESH_HPP

#include <cstdio>
#include "yafel_globals.hpp"
#include "mesh/GenericMesh.hpp"
#include "output/VTKObject.hpp"
#include "utils/ElementType.hpp"
#include "utils/ElementVtkType.hpp"

YAFEL_NAMESPACE_OPEN

template<typename MT, unsigned NSD>
class VTKMesh : public VTKObject {
  
private:
  using size_type = typename GenericMesh<MT,NSD>::size_type;
  const GenericMesh<MT,NSD> & M;



public:
  VTKMesh(const GenericMesh<MT,NSD> &m);
  void write(FILE *fp);
};



/*
 * Implementation
 */
template<typename MT, unsigned NSD>
VTKMesh<MT,NSD>::VTKMesh(const GenericMesh<MT,NSD> &m) :
  VTKObject(NONE, VTKMESH, std::string("Mesh")), M(m)
{}


template<typename MT, unsigned NSD>
void VTKMesh<MT,NSD>::write(FILE *fp) {
  
  
  //count cells to be written
  size_type nCells = 0;
  for(size_type i=0; i < M.n_elements(); ++i) {
    ElementType et = M.element_type(i);
    if(et != ElementType::NULL_ELEMENT)
      ++nCells;
  }
  
  //write header
  fprintf(fp, "%s\n", "<UnstructuredGrid>");
  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
	  M.n_nodes(), nCells );
  
  //write point data
  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  
  for(size_type nodeNum=0; nodeNum < M.n_nodes(); ++nodeNum) {
    fprintf(fp, "%f %f %f\n",
	    M.node(nodeNum)(0),
	    (NSD >= 2) ? M.node(nodeNum)(1) : 0.0,
	    (NSD >= 3) ? M.node(nodeNum)(2) : 0.0);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Points>\n");
  
  //write cell data
  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n");

  for(size_type elnum=0; elnum < M.n_elements(); ++elnum) {
    ElementType et = M.element_type(elnum);

    if(et==ElementType::NULL_ELEMENT)
      continue;
    
    switch(et) {
    case ElementType::LINEAR_LINE:
    case ElementType::QUADRATIC_LINE:
    case ElementType::CUBIC_LINE:
      fprintf(fp, "%d ", M.element(elnum)[0]);
      for (size_type i=2; i<M.element(elnum).size(); ++i) {
        fprintf(fp, "%d ", M.element(elnum)[i]);
      }
      fprintf(fp, "%d ", M.element(elnum)[1]);
      fprintf(fp, "\n");
      break;

    default:
      for(size_type i=0; i<M.element(elnum).size(); ++i) {
        fprintf(fp, "%d ", M.element(elnum)[i]);
      }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "</DataArray>\n");
  
  fprintf(fp, "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
  int offset = 0;
  for(size_type elnum=0; elnum<M.n_elements(); ++elnum) {
    ElementType et = M.element_type(elnum);
    if(et==ElementType::NULL_ELEMENT)
      continue;
    
    offset += M.element(elnum).size();
    fprintf(fp, "%d\n", offset);
  }
  fprintf(fp, "</DataArray>\n");
  
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(size_type elnum=0; elnum < M.n_elements(); ++elnum) {
    ElementType et = M.element_type(elnum);
    if(et == ElementType::NULL_ELEMENT)
      continue;
    
    fprintf(fp, "%d\n", ElementType_to_ElementVtkType(et));
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");
  
} // end write(FILE* fp)


YAFEL_NAMESPACE_CLOSE

#endif
