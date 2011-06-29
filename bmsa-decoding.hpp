/*
 * bmsa-decoding.hpp
 *
 *  Created on: 29.06.2011
 *      Author: ulysses
 */

#ifndef BMSA_DECODING_HPP_
#define BMSA_DECODING_HPP_

#include <array>
#include <vector>

#include "Point.hpp"
#include "bmsa.hpp"
#include "mv_poly.hpp"

namespace mv_poly {

template<typename T, int N>
struct CurvePoint {
    std::array<T, N> data;
};

template<
    int Dim,
    typename Field,
    template <typename PointImpl> class OrderPolicy = GradedAntilexMonomialOrder
>
class BMSDecoding :
    BMSAlgorithm<typename MVPolyType<Dim, Field>::type, OrderPolicy> {

public:

    typedef std::vector< CurvePoint<Field, Dim> > CurvePointsCollection;

    typedef std::vector<Field> FieldElemsCollection;

private:

    CurvePointsCollection points;

    // compute roots of elements in F
    CurvePointsCollection
    getErrorLocations() {

    }

    FieldElemsCollection
    getErrorValues(
            FieldElemsCollection const & r,
            CurvePointsCollection const & locations) {

    }

    void computeErrorLocatorPolynomials() {

    }

public:

    BMSDecoding(
            CurvePointsCollection const & points) : points(points) {}

    FieldElemsCollection decode(FieldElemsCollection const & r) {
        computeErrorLocatorPolynomials();
        CurvePointsCollection locations = getErrorLocations();
        FieldElemsCollection values = getErrorValues(r, locations);
        // correct r with locators and values...
    }
}; // BMSDecoding

} // namespace mv_poly

#endif /* BMSA_DECODING_HPP_ */
