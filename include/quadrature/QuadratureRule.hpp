#ifndef _YAFEL_QUADRATURERULE_HPP
#define _YAFEL_QUADRATURERULE_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \class QuadratureRule
 * \brief Class to represent a quadrature rule.
 *
 */
class QuadratureRule
{
public:
    enum class QuadratureType
    {
        GAUSS_LEGENDRE,
        GAUSS_LOBATTO
    };

    /**
     * Constructs an empty quadrature rule
     */
     QuadratureRule() : nodes(), weights() {}


    /**
     * constructs a 1D rule for \f$ \x \in [-1,1] \f$.
     * Higher-dimension quadrature rules should be constructed externally
     * and assigned into the "weights" and "nodes" members.
     * @param nPoints - number of integration points to use
     * @param quadratureType
     */
    QuadratureRule(int nPoints, QuadratureType quadratureType = QuadratureType::GAUSS_LEGENDRE);


    /**
     * Construct a quadrature rule for a tensor-product element
     * @param qt Quadrature Type
     * @param topoDim Dimension of tensor product rule {1,2,3}
     * @param polyOrder Polynomial order to integrate exactly (often useful to use 2x element polynomial interpolation)
     * @return
     */
    static QuadratureRule make_tensor_product(QuadratureType qt, int topoDim, int polyOrder);

    void get_triangle_quadrature(int porder);

    std::vector<coordinate<>> nodes;
    std::vector<double> weights;

private:
    void make_gauss_lobatto_1D(int nPoints);

    void make_gauss_legendre_1D(int nPoints);
};


YAFEL_NAMESPACE_CLOSE

#endif
