//
// Created by tyler on 9/12/17.
//

#ifndef YAFEL_CELL_CENTROIDS_HPP
#define YAFEL_CELL_CENTROIDS_HPP

#include "yafel.hpp"

YAFEL_NAMESPACE_OPEN

template<int NSD>
auto CellCentroidsVolumes(const DoFManager &dofm)
{


    TaskScheduler TS(4);

    std::vector<ElementFactory> EFs(TS.workers.size());
    //ElementFactory EF;
    std::vector<coordinate<>> centroids(dofm.nCells(), coordinate<>(0));
    std::vector<double> cell_measure(dofm.nCells(), 0);

    std::mutex mtx_1;

    struct thread_local_variables {
        ElementFactory EF;
    };


    std::vector<thread_local_variables> local_variables(TS.workers.size(), {ElementFactory(1)});

    auto loop_body = [&local_variables, &centroids, &cell_measure, &mtx_1, &dofm](std::size_t elnum) {

        auto &EF = local_variables[worker_global::worker_id].EF;
        double Volume{0};
        coordinate<> centroid_e(0);

        auto et = dofm.element_types[elnum];
        if (et.topoDim != NSD) {
            return;
        }
        //Element &E = EFs[worker_global::worker_id].getElement(et);
        auto &E = EF.getElement(et);

        for (int qpi = 0; qpi < E.nQP(); ++qpi) {
            E.template update<NSD>(elnum, qpi, dofm);

            //get quadrature point coordinate
            coordinate<> xqp(0);
            for (auto A = 0; A < E.shapeValues[qpi].rows(); ++A) {
                xqp += dofm.dof_nodes[E.globalNodes[A]] * E.shapeValues[qpi](A);
            }

            centroid_e += xqp * E.jxw;
            Volume += E.jxw;
        }

        centroids[elnum] = centroid_e / Volume;
        cell_measure[elnum] = Volume;

    };

    auto futures = parfor(0, centroids.size(), loop_body, TS, 1);
    for (auto &f : futures) {
        f.wait();
    }


    return std::make_pair(std::move(centroids), std::move(cell_measure));
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_CELL_CENTROIDS_HPP
