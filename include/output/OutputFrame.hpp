//
// Created by tyler on 4/11/17.
//

#ifndef YAFEL_OUTPUTFRAME_HPP
#define YAFEL_OUTPUTFRAME_HPP

#include "yafel_globals.hpp"
#include "output/OutputData.hpp"
#include "output/OutputBackend.hpp"
#include <memory>


YAFEL_NAMESPACE_OPEN

class OutputFrame
{
public:
    inline OutputFrame(OutputMesh &om) : outputMesh(&om) {}

    OutputMesh *outputMesh;
    std::vector<std::shared_ptr<OutputData>> point_data;
    std::vector<std::shared_ptr<OutputData>> cell_data;
    double time{0};


    inline void addData(std::shared_ptr<OutputData> data)
    {
        if (data->dataLocation == OutputData::DataLocation::Point)
            point_data.push_back(data);
        if (data->dataLocation == OutputData::DataLocation::Cell)
            cell_data.push_back(data);
    }
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTFRAME_HPP
