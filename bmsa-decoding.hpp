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

#include <glog/logging.h>

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
    typename ECCodeParams/*,
    template <typename> class OrderPolicy = GradedAntilexMonomialOrder*/
>
class BMSDecoding {

public:

    typedef typename ECCodeParams::BasisElem BasisElem;

    typedef typename ECCodeParams::CurvePoint CurvePoint;

    typedef std::vector< CurvePoint > CurvePointsCollection;

    typedef std::vector<Field> FieldElemsCollection;

private:

    /******************** Private typedef's *********************/

    typedef std::map< BasisElem, Field> SyndromeType;

    typedef typename ECCodeParams::OrderPolicyHolder OrderPolicyHolder;

    typedef BMSAlgorithm< SyndromeType,
            typename MVPolyType<Dim, Field>::type,
            OrderPolicyHolder::template impl
        > BmsaT;

    typedef typename BmsaT::PolynomialCollection PolynomialCollection;

    typedef std::vector<int> ErrorPositions;

    /******************** Private fields *********************/
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

    PolynomialCollection
    computeErrorLocatorPolynomials(FieldElemsCollection const & received) {
        SyndromeType syn;
        auto basis = ECCodeParams::getCodeBasis(l);

        // **********  logging
        ostringstream log_oss;
        std::copy(basis.begin(), basis.end(),
                std::ostream_iterator<BasisElem>(log_oss, " "));
        LOG(INFO) << log_oss.str();
        // **********  ENF OF logging

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

        // **********  logging
        log_oss.str("");
        std::for_each(syn.begin(), syn.end(),
                [&log_oss](typename SyndromeType::value_type const & s) {
            log_oss << makeNtlPowerPrinter(s.second) << " ";

        });
        LOG(INFO) << log_oss.str();
        // **********  ENF OF logging

        BmsaT bmsa(syn, ++basis.back());
        auto minset = bmsa.computeMinimalSet();

        // **********  logging
        //log_oss.str("");
        for_each(minset.begin(), minset.end(),
                [/*&log_oss*/](typename BmsaT::PolynomialCollection::value_type const & p) {
                    LOG(INFO) <<
                            makePowerPrinter< OrderPolicyHolder::template impl >(p)
                            << std::endl;
                }
        );
//        LOG(INFO) << log_oss.str();
        // **********  ENF OF logging

        return minset;
//        while(true) { // TODO: define stop condition
//            frmv();
//            ++k; // TODO: ++ should be correctly implemented for given OrderPolicy
//        }
    }

public:

    BMSDecoding(
            size_t l)
    : l(l) {
        curvePoints = ECCodeParams::getRationalPoints();
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
