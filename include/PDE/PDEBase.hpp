//
// Created by tyler on 5/23/17.
//

#ifndef YAFEL_PDEBASE_HPP
#define YAFEL_PDEBASE_HPP

/**
 * \class PDEBase
 * \brief Base class for all classes implementing PDE-specific physics
 */
template<int NSD>
struct PDEBase {
    static constexpr int nsd() { return NSD; }

    /**
     * \class PDEDefaultData
     * \brief Default class of PDEData, which will hold
     * data used by specific physics models (eg material models),
     * and will be instantiated, one per quadrature point
     * by the FESystem.
     */
    struct PDEDefaultData {};


    using PDEData = PDEDefaultData;
};

#endif //YAFEL_PDEBASE_HPP
