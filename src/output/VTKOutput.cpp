#include "output/VTKOutput.hpp"
#include <cstdlib>
#include <cstdio>
#include <stdexcept>

YAFEL_NAMESPACE_OPEN

VTKOutput::VTKOutput() :
    vtkmesh(nullptr),
    pointData(),
    cellData()
{}

void VTKOutput::addVTKObject(VTKObject *VO) {
    switch(VO->getObjectType()) {
    case VTKObject::VTKMESH:
        vtkmesh = VO;
        break;
    case VTKObject::VTKPOINTDATA:
        pointData.push_back(VO);
        break;
    case VTKObject::VTKCELLDATA:
        cellData.push_back(VO);
        break;
    }
}

void VTKOutput::clearData() {
    pointData.clear();
    cellData.clear();
}

void VTKOutput::clearMesh() {
    vtkmesh = NULL;
}

void VTKOutput::write(const std::string &fname) {
  
    FILE *fp = fopen(fname.c_str(), "w");
    if(fp == NULL) {
        std::runtime_error("VTKOutput::write() : Could not open file.");
    }

    if(vtkmesh == NULL) {
        std::runtime_error("VTKOutput::write() : no mesh assigned.");
    }
  
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
    vtkmesh->write(fp);
  
    //write point data
    if(pointData.size() > 0) {
        fprintf(fp, "<PointData>\n");
        for(unsigned i=0; i<pointData.size(); ++i) {
            pointData[i]->write(fp);
        }
        fprintf(fp, "</PointData>\n");
    }

    //write cell data
    if(cellData.size() > 0) {
        fprintf(fp, "<CellData>\n");
        for(unsigned i=0; i<cellData.size(); ++i) {
            cellData[i]->write(fp);
        }
        fprintf(fp, "</CellData>\n");
    }

    //conclude file
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>");

    fclose(fp);
  
    fprintf(stdout, "VTK Output successfully written to %s\n", fname.c_str());
}

YAFEL_NAMESPACE_CLOSE
