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

#include <boost/iterator/transform_iterator.hpp>

#include "Point.hpp"
#include "bmsa.hpp"
#include "mv_poly.hpp"
#include "CurveArithmetic.hpp"

namespace mv_poly {

//template<typename T, int N>
//struct CurvePoint {
//    std::array<T, N> data;
//};

template<
    int Dim,
    typename Field,
    typename ECCodeTraits/*,
    template <typename> class OrderPolicy = GradedAntilexMonomialOrder*/
>
class BMSDecoding {

public:

    typedef typename ECCodeTraits::BasisElem BasisElem;

    typedef typename ECCodeTraits::CurvePoint CurvePoint;

    typedef std::vector< CurvePoint > CurvePointsCollection;

    typedef std::vector<Field> FieldElemsCollection;

private:

//    typedef Point<Dim, OrderPolicy> LatticePoint;

    typedef std::map< BasisElem, Field> SyndromeType;

    size_t l;

    CurvePointsCollection curvePoints;



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

    void computeErrorLocatorPolynomials(FieldElemsCollection const & received) {
        SyndromeType syn;
        auto basis = ECCodeTraits::getCodeBasis(l);

        using namespace std::tr1::placeholders;
        auto syndromComponentAtBasisElem =
            [this,&received](BasisElem const & be) -> typename SyndromeType::value_type {
                auto tit = boost::make_transform_iterator(this->curvePoints.begin(),
                        std::tr1::bind(
                            computeMonomAtPoint<Field, BasisElem, CurvePoint>,
                            be,
                            _1));
                return typename SyndromeType::value_type(
                        be,
                        std::inner_product(received.begin(), received.end(),
                            tit, FieldElemTraits<Field>::addId()));

        };
        std::transform(basis.begin(), basis.end(),
                std::inserter(syn, syn.begin()), syndromComponentAtBasisElem);
        std::for_each(syn.begin(), syn.end(),
                [](typename SyndromeType::value_type const & s) {
            std::cout << makeNtlPowerPrinter(s.second) << " ";

        });
        //        auto sn = f(basis[0]);
        //        std::cout << "syn[0]: " << //makeNtlPowerPrinter(
        //                   sn/*)*/ << std::endl;

        //BMSAlgorithm< SyndromeType, MVPolyType<Dim, Field>::type >
            // bmsa(syn, synLength);

//        while(true) { // TODO: define stop condition
//            frmv();
//            ++k; // TODO: ++ should be correctly implemented for given OrderPolicy
//        }
    }

public:

    BMSDecoding(
            size_t l)
    : l(l) {
        curvePoints = ECCodeTraits::getRationalPoints();
    }

    FieldElemsCollection decode(FieldElemsCollection const & r) {
        FieldElemsCollection result;
        computeErrorLocatorPolynomials(r);
//        CurvePointsCollection locations = getErrorLocations();
//        FieldElemsCollection values = getErrorValues(locations);
        // correct r with locators and values...
        return result;
    }
}; // BMSDecoding

} // namespace mv_poly

#endif /* BMSA_DECODING_HPP_ */
