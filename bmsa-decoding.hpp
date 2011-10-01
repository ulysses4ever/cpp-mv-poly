/*
 * bmsa-decoding.hpp
 *
 *  Created on: 29.06.2011
 *      Author: ulysses
 */

#ifndef BMSA_DECODING_HPP_
#define BMSA_DECODING_HPP_

#include <array>
#include <map>
#include <vector>

#include "Point.hpp"
#include "bmsa.hpp"
#include "mv_poly.hpp"
#include "CurveArithmetic.hpp"

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
class BMSDecoding {

public:

    typedef std::vector< CurvePoint<Field, Dim> > CurvePointsCollection;

    typedef std::vector<Field> FieldElemsCollection;

private:

    typedef std::map< Point<Dim, OrderPolicy>, Field> SyndromeType;

    size_t a;

    CurvePointsCollection points;

    //FieldElemsCollection * r = 0;

    //BMSAlgorithm<typename MVPolyType<Dim, Field>::type, OrderPolicy> * bmsa = 0;

    // compute roots of elements in F
    CurvePointsCollection
    getErrorLocations() {

    }

    FieldElemsCollection
    getErrorValues(
            FieldElemsCollection const & r,
            CurvePointsCollection const & locations) {

    }

    // Feng-Rao majority voting
    void frmv() {

    }

    void computeErrorLocatorPolynomials(FieldElemsCollection const & r) {
        int i = 0;
        Point<Dim, OrderPolicy> k, synLength;
        SyndromeType syn;
        // compute known syndrome length on the basis of code parameter 'a'
        for (int i = 0; i < a; ++i) {
            syn[k] = in_prod(expPoints(k), r);
            ++k;
        }
        //BMSAlgorithm< SyndromeType, MVPolyType<Dim, Field>::type >
            // bmsa(syn, synLength);

        while(true) { // TODO: define stop condition
            frmv();
            ++k; // TODO: ++ should be correctly implemented for given OrderPolicy
        }
    }

public:

    BMSDecoding(
            size_t a,
            CurvePointsCollection const & points)
    : a(a), points(points) {}

    FieldElemsCollection decode(FieldElemsCollection const & r) {
        computeErrorLocatorPolynomials(r);
        CurvePointsCollection locations = getErrorLocations();
        FieldElemsCollection values = getErrorValues(locations);
        // correct r with locators and values...
    }
}; // BMSDecoding

} // namespace mv_poly

#endif /* BMSA_DECODING_HPP_ */
