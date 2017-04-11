//
// Created by tyler on 4/11/17.
//

#ifndef YAFEL_OUTPUTFRAME_HPP
#define YAFEL_OUTPUTFRAME_HPP

#include "yafel_globals.hpp"
#include "output/OutputData.hpp"
#include "output/OutputBackend.hpp"

YAFEL_NAMESPACE_OPEN

class OutputFrame
{
public:
    const OutputMesh *outputMesh;
    std::vector<OutputData*> point_data;
    std::vector<OutputData*> cell_data;
    double time;

    //void addPointData();
    //void addCellData();

};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTFRAME_HPP
