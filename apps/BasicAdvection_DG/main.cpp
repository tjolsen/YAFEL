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
    std::vector<Eigen::PartialPivLU<Eigen::MatrixXd>> inverse_mass_matrices;
    bool mass_constructed{false};

    static constexpr int nsd()
    { return NSD; }

    static double initial_condition(const coordinate<> &x)
    {

        double r0 = 0.15;
        double r = norm(x - coordinate<>{1.5, 1});
        return std::exp(-(r / r0) * (r / r0));
    }

    static Tensor<NSD, 1> convection_velocity(coordinate<> x, double)
    {

        //double omega = 0.10;
        //coordinate<> dx = x - coordinate<>{1, 1, x(2)};

        //double r = norm(dx);
        //double theta = std::atan2(dx(1), dx(0));

        //Tensor<NSD, 1> v{-r * omega * std::sin(theta), r * omega * std::cos(theta)};
        //return v;
        return {-1, .5};
    };

    static double source(const coordinate<> &, double)
    {
        return 0.0;
    }

    static void LocalResidual(const Element &E, int qpi, double time,
                              Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &U_el,
                              Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> &R_el)
    {
        coordinate<> xqp(0);
        for (int i = 0; i < E.nQP(); ++i) {
            xqp += E.shapeValues[qpi](i) * E.localMesh.getGeometryNodes()[i];
        }
        auto c = convection_velocity(xqp, time);

        double U{0};
        for(int A=0; A<U_el.rows(); ++A) {
            U += U_el(A)*E.shapeValues[qpi](A);
        }

        for (int i = 0; i < R_el.rows(); ++i) {
            auto gradW = make_TensorMap<NSD, 1>(&E.shapeGrad(i, 0));
            auto val = (U * dot(gradW, c)
                        + source(xqp, time) * E.shapeValues[qpi](i));
            R_el(i) += val * E.jxw;
        }
    }

    template<typename T>
    static void LocalMass(const Element &E, int qpi, double,
                          Eigen::MatrixBase<T> &M_el)
    {
        M_el += (E.jxw * E.shapeValues[qpi]) * E.shapeValues[qpi].transpose();
    }


    static double BoundaryFlux(Tensor<NSD, 1> n, coordinate<> x, double t, double U)
    {
        auto vdotn = dot(n, convection_velocity(x, t));
        if(vdotn > 0) {
            return vdotn*U;
        }
        else {
            return 0;
        }
    }

    static double Flux(Tensor<NSD, 1> n, coordinate<> x, double t, double Uplus, double Uminus)
    {
        auto vdotn = dot(n, convection_velocity(x, t));
        return 0.5 * (vdotn * (Uplus + Uminus) + std::abs(vdotn) * (Uplus - Uminus));
    }

};


int main()
{
    constexpr int NSD = 2;
    Mesh M("mesh.msh");
    M.buildInternalFaces();
    int p = 2;
    int dofpn = 1;

    AdvectionPhysics<NSD> AP;

    double Tfinal = 1;
    double dt = .001;
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

        DGAssembly<AdvectionPhysics<NSD>>(feSystem, AP,
                                          {AssemblyRequirement::Residual,
                                           AssemblyRequirement::DtMass});
        feSystem.getSolution() += dt * feSystem.getGlobalResidual();

        if (ti % output_freq == 0) {
            simulationOutput.captureFrame(feSystem);
        }
    }

}