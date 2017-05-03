//
// Created by tyler on 5/2/17.
//

#include "element/ElementFactory.hpp"
#include "assembly/DGAssembly.hpp"
#include "fe_system/FESystem.hpp"
#include "output/SimulationOutput.hpp"

using namespace yafel;

template<int NSD>
struct AdvectionPhysics
{
    static constexpr int nsd() { return NSD; }

    static double initial_condition(const coordinate<> &x)
    {
        /*
        double r0 = 0.15;
        double r = norm(x - coordinate<>{1.5, 1});
        return std::exp(-(r / r0) * (r / r0)) * (r <= .25);
         */

        return std::exp(-std::abs(x(0)-1));
    }

    static Tensor<NSD, 1> convection_velocity(coordinate<> x, double)
    {

        //double omega = 1.0;
        //coordinate<> dx = x - coordinate<>{1, 1, x(2)};

        //double r = norm(dx);
        //double theta = std::atan2(dx(1), dx(0));

        //Tensor<NSD, 1> v{-r * omega * std::sin(theta), r * omega * std::cos(theta)};
        return {1,0};
    };

    static double source(const coordinate<> &, double)
    {
        return 0.0;
    }

    static void LocalResidual(const Element &E, int qpi, double time,
                              Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &U_el,
                              Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {
        coordinate<> xqp;
        for (int i = 0; i < E.nQP(); ++i) {
            xqp = xqp + E.shapeValues[qpi](i) * E.localMesh.getGeometryNodes()[i];
        }
        auto c = convection_velocity(xqp, time);


        for (int i = 0; i < R_el.rows(); ++i) {
            R_el(i) += (U_el(i) * dot(make_TensorMap<NSD, 1>(&E.shapeGrad(i, 0)), c)
                        + source(xqp, time) * E.shapeValues[qpi](i)) * E.jxw;
        }
    }

    static void LocalMass(const Element &E, int qpi, double,
                          Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &M_el)
    {
        M_el += (E.jxw * E.shapeValues[qpi]) * E.shapeValues[qpi].transpose();
    }


    static double Flux(Tensor<NSD, 1> n, coordinate<> x, double t, double Uplus, double Uminus)
    {
        auto vdotn = dot(n, convection_velocity(x, t));
        return 0.5 * (vdotn * (Uplus+Uminus) + std::abs(vdotn)*(Uplus-Uminus));
    }

};


int main()
{
    constexpr int NSD = 2;
    Mesh M("mesh.msh");
    M.buildInternalFaces();
    int p = 3;
    int dofpn = 1;


    double Tfinal = .1;
    double dt = .0001;
    int output_freq = 1;


    DoFManager dofm(M, DoFManager::ManagerType::DG, p, dofpn);

    FESystem feSystem(dofm);

    auto &U = feSystem.getSolution();

    //set initial condition
    for (int i = 0; i < dofm.dof_nodes.size(); ++i) {
        U(i) = AdvectionPhysics<NSD>::initial_condition(dofm.dof_nodes[i]);
    }

    SimulationOutput simulationOutput("output", BackendType::HDF5);

    simulationOutput.captureFrame(feSystem);

    int NT = static_cast<int>(Tfinal / dt);
    for (int ti = 1; ti <= NT; ++ti) {
        std::cout << ti << " / " << NT << std::endl;
        feSystem.currentTime() = ti * dt;

        DGAssembly<AdvectionPhysics<NSD>>(feSystem,
                                          {AssemblyRequirement::Residual,
                                           AssemblyRequirement::DtMass});
        feSystem.getSolution() += dt * feSystem.getGlobalResidual();

        if (ti % output_freq == 0) {
            simulationOutput.captureFrame(feSystem);
        }
    }

}