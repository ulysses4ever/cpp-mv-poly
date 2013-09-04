/*
 * bmsa-decoding.hpp
 *
 *  Created on: 29.06.2011
 *      Author: ulysses
 */

#ifndef BMSA_DECODING_HPP_
#define BMSA_DECODING_HPP_

#include <array>
#include <functional>
#include <map>
#include <vector>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/iterator_range.hpp>
//#include <boost/bind.hpp>

#include <glog/logging.h>

#include "Point.hpp"
#include "bmsa.hpp"
#include "mv_poly.hpp"
#include "CurveArithmetic.hpp"

namespace mv_poly {

template<
    int Dim,
    typename ECCodeParams
>
class BMSDecoding {

public:

    typedef typename ECCodeParams::Field Field;

    typedef typename ECCodeParams::BasisElem BasisElem;

    typedef BasisElem PointT;

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

    // compute common roots of elements in F
    ErrorPositions
    getErrorLocations(PolynomialCollection const & polys) {
        typedef typename PolynomialCollection::value_type PolyT;
        typedef CoefficientTraits<typename PolyT::CoefT> CoefTr;
        ErrorPositions result;
        int idx = 0;
        std::for_each(curvePoints.begin(), curvePoints.end(),
                [&result,&idx,&polys](CurvePoint const & cp) {
                    bool root_for_all = std::accumulate(polys.begin(),
                            polys.end(), true,
                            [&cp](bool val, PolyT const & p){
                                return val && p(cp) == CoefTr::addId();
                            });
                    if (root_for_all)
                        result.push_back(idx);
                    ++idx;

        });

        // **********  logging
        std::ostringstream log_oss;
        std::copy(result.begin(), result.end(),
                std::ostream_iterator<int>(log_oss, " "));
        LOG(INFO) << "error positions: " << log_oss.str();
        // **********  ENF OF logging

        return result;
    }

    FieldElemsCollection
    getErrorValues(
            FieldElemsCollection const & r,
            CurvePointsCollection const & locations) {

    }

    // Feng-Rao majority voting
    void frmv(BmsaT & bmsa) {
        typedef std::map<PointT, PointT> PointToPointMap;
        PointToPointMap voters; // Г_k with support points from \sigma_k attached
        PointT k = bmsa.getSeqLen();     // k - kurrent step :)

        // 1. Compute Г_k
        for (
             auto it = bmsa.getF().begin();
             it != bmsa.getF().end();
             ++it
             ) {
            PointT const & s = it->first;
            for (PointT t = s; t < k; ++t) {
                // we assume tha going from s to k following monomial order
                // yields all points t s.t. s <= t <= k where <= stands for
                // by-coordinate comparison (natural partial order)
                bool l1 = byCoordinateLess(s, t);
                bool l2 = byCoordinateLess(t, k);
                bool g = byCoordinateGreaterThenAny(k - t, bmsa.getF() | boost::adaptors::map_keys);
                if (
                        l1 &&
                        l2 &&
                        g) {
                    voters[t] = s;
                }
            }
        }
        // **********  logging (Г_k)
        LOG(INFO) << "Gamma_k has " << voters.size() << " elements";
        std::ostringstream log_oss;
        std::copy(voters.begin(), voters.end(),
                std::ostream_iterator<typename PointToPointMap::value_type>(log_oss, " "));
        LOG(INFO) << "Gamma_k: " << log_oss.str();
        // **********  ENF OF logging


        // 2. Compute votes and choose winner

    }

    PolynomialCollection
    computeErrorLocatorPolynomials(FieldElemsCollection const & received) {
        SyndromeType syn;
        auto basis = ECCodeParams::getCodeBasis(l);

        // **********  logging
        std::ostringstream log_oss;
        std::copy(basis.begin(), basis.end(),
                std::ostream_iterator<BasisElem>(log_oss, " "));
        LOG(INFO) << "Basis for evaluation code: " << log_oss.str();
        // **********  ENF OF logging

        // computing "known" syndroms
        //using namespace boost;//
        using namespace std::placeholders;
        auto syndromComponentAtBasisElem =
            [this,&received](BasisElem const & be) -> typename SyndromeType::value_type {
                auto tit = boost::make_transform_iterator(this->curvePoints.begin(),
                        std::bind(
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
        // END OF computing "known" syndroms

        // **********  logging
        log_oss.str("");
        std::for_each(syn.begin(), syn.end(),
                [&log_oss](typename SyndromeType::value_type const & s) {
                    log_oss << makeNtlPowerPrinter(s.second) << " ";
                }
        );
        LOG(INFO) << "Known syndroms: " << log_oss.str();
        // **********  ENF OF logging

        // compute error locators for the known syndroms
        BmsaT bmsa(syn, ++basis.back());
        auto minset = bmsa.computeMinimalSet();

        // **********  logging
        //log_oss.str("");
        for_each(minset.begin(), minset.end(),
                [/*&log_oss*/](typename BmsaT::PolynomialCollection::value_type const & p) {
                    LOG(INFO) << "Error locator: " <<
                            makePowerPrinter< OrderPolicyHolder::template impl >(p)
                            << std::endl;
                }
        );
//        LOG(INFO) << log_oss.str();
        // **********  ENF OF logging

//        std::cout << "Hi!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
//        while(false) { // TODO: define stop condition
            frmv(bmsa);
//            ++bmsa; // TODO: ++ should be correctly implemented for given OrderPolicy
//        }
        return minset;
    }

public:

    BMSDecoding(
            size_t l)
    : l(l) {
        curvePoints = ECCodeParams::getRationalPoints();
    }

    FieldElemsCollection decode(FieldElemsCollection const & r) {
        FieldElemsCollection result;
        ErrorPositions locations = getErrorLocations(
                computeErrorLocatorPolynomials(r));
//        FieldElemsCollection values = getErrorValues(locations);
        // correct r with locators and values...
        return result;
    }
}; // BMSDecoding

} // namespace mv_poly

#endif /* BMSA_DECODING_HPP_ */
