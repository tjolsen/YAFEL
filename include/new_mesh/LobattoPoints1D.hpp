//
// Created by tyler on 3/14/17.
//

#ifndef YAFEL_LOBATTOPOINTS1D_HPP
#define YAFEL_LOBATTOPOINTS1D_HPP

#include "yafel_globals.hpp"
#include <vector>


YAFEL_NAMESPACE_OPEN

/**
 * \brief Make the Lobatto Interpolation points on the biunit line \f$ x\in [-1,1] \f$
 *
 * Returns a vector of polynomialOrder+1 points which will include the enpoints
 *
 * @param polynomialOrder Order of interpolating polynomials
 * @return vector of points on \f$ x\in[-1, 1] \f$
 */
std::vector<double> make_LobattoPoints_1D(int polynomialOrder);

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_LOBATTOPOINTS1D_HPP
