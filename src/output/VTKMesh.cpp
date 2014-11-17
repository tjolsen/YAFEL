#include "output/VTKMesh.hpp"

YAFEL_NAMESPACE_OPEN

VTKMesh::VTKMesh(ElementFactory &EF) : 
  VTKObject(NONE, VTKMESH, std::string("Mesh")), EFp(&EF) {}

void VTKMesh::write(FILE *fp) {
  //get pointer to Mesh
  Mesh *Mp = EFp->getMesh();

  int nCells = 0;
  for(unsigned i=0; i < Mp->get_n_elems(); ++i) {
    Element *e = EFp->getElement(i);
    if(e != NULL)
      nCells += 1;
  }
  
  fprintf(fp, "%s\n", "<UnstructuredGrid>");
  fprintf(fp, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
	  Mp->get_n_nodes(), nCells );
  
  //write point data
  fprintf(fp, "<Points>\n");
  fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  
  for(unsigned nodeNum=0; nodeNum < Mp->get_n_nodes(); ++nodeNum) {
    fprintf(fp, "%f %f %f\n",
	    Mp->nodal_coords[nodeNum](0),
	    Mp->nodal_coords[nodeNum](1),
	    Mp->nodal_coords[nodeNum](2));
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Points>\n");
  
  //write cell data
  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n");

  for(unsigned elnum=0; elnum < Mp->get_n_elems(); ++elnum) {
    Element *e = EFp->getElement(elnum);
    if(e==NULL)
      continue;

    for(unsigned i=0; i<Mp->elements[elnum].size(); ++i) {
      fprintf(fp, "%d ", Mp->elements[elnum][i]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  
  fprintf(fp, "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
  int offset = 0;
  for(unsigned elnum=0; elnum<Mp->get_n_elems();++elnum) {
    Element *e = EFp->getElement(elnum);
    if(e==NULL)
      continue;
    
    offset += e->nodes_per_el;
    fprintf(fp, "%d\n", offset);
  }
  fprintf(fp, "</DataArray>\n");

  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(unsigned elnum=0; elnum < Mp->get_n_elems(); ++elnum) {
    Element *e = EFp->getElement(elnum);
    if(e == NULL)
      continue;

    fprintf(fp, "%d\n", e->vtk_type);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");
}

YAFEL_NAMESPACE_CLOSE
