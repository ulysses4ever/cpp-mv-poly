/**
 * @file bmsa.hpp
 *
 * Implementation of Berlekamp—Massey—Sakata algorithm.
 *
 *  @date 2010-07-27
 *  @author Artem Pelenitsyn
 */

#ifndef BMSA_HPP_
#define BMSA_HPP_

#include <iostream>
#include <iomanip>
#include <iterator>
#include <list>
#include <map>
#include <utility>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <glog/logging.h>

#include <tr1/functional>

#include "mv_poly.hpp"
#include "Utilities.hpp"
#include "Point.hpp"
#include "CoefficientTraits.hpp"

namespace mv_poly {

template<
    typename SeqT,
    typename PolynomialT = SeqT,
    template <typename PointImpl> class OrderPolicy = GradedAntilexMonomialOrder
>
class BMSAlgorithm {
public:
    static const int Dim = PolynomialT::VAR_CNT;

    typedef Point<Dim, OrderPolicy> PointT;

    typedef std::map< PointT, PolynomialT > PointPolyMap;

    typedef std::list< PolynomialT > PolynomialCollection;

    typedef std::list< PointT > PointCollection;

private:

    typedef typename PolynomialT::CoefT CoefT;

    typedef std::map< PointT, CoefT > PointCoefMap;

    PointPolyMap G;

    const CoefT ZERO;

    PointT seqLen;

    PointCollection oldDeltaPoints;

    PointPolyMap F;

    SeqT & seq;

public:

    void infoUpdate(PointT const & k) {
        using std::tr1::bind;
        using std::tr1::cref;
        using namespace std::tr1::placeholders;
        using std::cout;
        using std::endl;
        LOG(INFO) << "Enter InfoUpdate(k), k = " << k << endl;

        PointCollection deltaPoints, sigmaPoints;
        PointCoefMap discr; // discrepancies

        //searching for candidates to form new deltaPoints
        LOG(INFO) << "traversing F (building new delta-set)" << endl;
        for (typename PointPolyMap::iterator fIt = F.begin();
                fIt != F.end(); ++fIt) {
            Point<Dim, OrderPolicy> const & degF = fIt->first;
            PolynomialT & f = fIt->second;
            if (byCoordinateLess(degF, k)) {
                CoefT b = conv(f, seq, degF, k);
                discr[degF] = b;
                Point<Dim, OrderPolicy> c = k - degF;
                LOG_IF(INFO, b != ZERO) << "\t d != 0, fallen:" << f;
                LOG_IF(INFO, b != ZERO) << "\t\tspan (c): " << c << endl;
                if  (b != ZERO &&
                        !byCoordinateLessThenAny(
                                c,
                                make_choose_point_iterator(G.begin()),
                                make_choose_point_iterator(G.end()))) {
                    deltaPoints.push_back(c);
                    LOG(INFO) << "\t\tc is a new delta-point";
                }
            }
        }
        
        deltaPoints.splice(deltaPoints.end(), oldDeltaPoints); // push_back all olds
        deltaPoints = getPartialMaximums(deltaPoints);
        sigmaPoints = getConjugatePointCollection(deltaPoints);

//            cout << "Delta-points:" << endl;
//            copy(deltaPoints.begin(), deltaPoints.end(),
//                    std::ostream_iterator< Point<Dim> >(cout, "\n"));
//            cout << "Sigma-points:" << endl;
//            copy(sigmaPoints.begin(), sigmaPoints.end(),
//                    std::ostream_iterator< Point<Dim> >(cout, "\n"));

        PointPolyMap newG;
        // forming new G
        LOG(INFO) << "forming new G (traversing new delta-set)" << endl;
        for (typename PointCollection::const_iterator
                cIt = deltaPoints.begin();
                cIt != deltaPoints.end(); ++cIt) {
            LOG(INFO) << "\tdelta-point: "  << *cIt << endl;
            typename PointPolyMap::const_iterator tmpIt = G.find(*cIt);
            if  (tmpIt != G.end()) {
                LOG(INFO) << "\t\told G can handle it: "  << *tmpIt << endl;
                newG.insert(*tmpIt);
            }   else {
                Point<Dim, OrderPolicy> s = k - *cIt;
                newG[*cIt] =
                        CoefficientTraits<CoefT>::multInverse(discr[s])
                        * F[s];
                LOG(INFO) << "\t\told f goes to new G: "  << newG[*cIt] << endl;
            }
        }

        PointPolyMap newF;
        // forming new F
        LOG(INFO) << "forming new F (traversing new sigma-points)" << endl;
        for (typename PointCollection::const_iterator
                tIt = sigmaPoints.begin();
                tIt != sigmaPoints.end(); ++tIt) {
            Point<Dim, OrderPolicy> const & t = *tIt;
            LOG(INFO) << "\tt: " << t << " ";
            
            // next find_if will find something for shure — I GUARANTEE IT
            auto fIt = std::find_if(
                    make_choose_point_iterator(F.begin()),
                    make_choose_point_iterator(F.end()),
                    bind(&byCoordinateLess<Dim, OrderPolicy>, _1, cref(t)));
            CHECK(fIt != make_choose_point_iterator(F.end()))
                << "CRITICAL: supporting f absent for given t: " << t;
                
            Point<Dim, OrderPolicy> s = *(fIt);
            Point<Dim, OrderPolicy> u = t - s; // valid Point as byCoordinateLess(s, t) == true
            CHECK(byCoordinateLess(s, t))
                << "CRITICAL: supporting f is not such that byCoordinateLess(deg(f), t), f: " 
                << *(fIt);

            LOG(INFO)
                << "\ts: " << s << " "
                << "u: " << u
                << endl;
            //LOG(INFO) << "F[s] << u: " << (F[s] << u) << endl;
            bool notJustIncreaseDegree = true; // some tricky flag to avoid goto
            if ((notJustIncreaseDegree = byCoordinateLess(t, k))) {
                // yes, I mean assignment in if condition
                const typename PointPolyMap::const_iterator cIt = (std::find_if(
                        make_choose_point_iterator(G.begin()),
                        make_choose_point_iterator(G.end()),
                        bind(&byCoordinateLess<Dim, OrderPolicy>, cref(k - t), _1))
                        ).base();
                if ((notJustIncreaseDegree = (cIt != G.end()))) {
                    // yes, I mean assignment at the top of if condition
                    Point<Dim, OrderPolicy> const & c = cIt->first;
                    newF[t] = (F[s] << u) - // Berlekamp formula
                        (discr[s] * G[c] << (c - (k - t)));
                    LOG(INFO) << "\tnew f via Berlekamp formula (deg is const)";
                }
            }
            if (!notJustIncreaseDegree) {
                newF[t] = F[s] << u;
                LOG(INFO) << "\tnew f via deg bump";
            }
        } // end of forming new F
        F = newF;
        G = newG;
        oldDeltaPoints.clear();
        oldDeltaPoints.splice(oldDeltaPoints.end(), deltaPoints);
            // this splice means: move contents of deltaPoints to oldDeltaPoints

        LOG(INFO) << "F: " << endl;
        copy(F.begin(), F.end(),
                std::ostream_iterator<typename PointPolyMap::value_type>(LOG(INFO), "\n\t"));
        LOG(INFO) << "G: " << endl << "\t";
        copy(G.begin(), G.end(),
                std::ostream_iterator<typename PointPolyMap::value_type>(LOG(INFO), "\n"));
        LOG(INFO) << endl;
    }

    BMSAlgorithm(
            SeqT & seq_,
            PointT const & seqLen_) :
                ZERO(CoefficientTraits<CoefT>::addId()),
                seqLen(seqLen_), seq(seq_) {
        F.insert(std::make_pair(PointT(), PolynomialT::getId()));
    }

    PolynomialCollection computeMinimalSet() {
        LOG(INFO) << "bmsa<" << Dim << ">::computeMinimalSet";
        // can't print seq as we SeqT can be various things 
         LOG(INFO) << mapToStr(seq);
                
        oldDeltaPoints.clear();
        // scanning input sequense step-by-step, following monomial order
        for (Point<Dim, OrderPolicy> k; k < seqLen; ++k) {
            infoUpdate(k);
        }
        return getPolynomialList();

    }

    PointPolyMap const & getF() const {
        return F;
    }

    PointCollection const & getDeltaPoints() const {
        return oldDeltaPoints;
    }

    PointT & getSeqLen() {
        return seqLen;
    }

    SeqT & getSeq() {
        return seq;
    }

private:
    PolynomialCollection getPolynomialList() {
        PolynomialCollection result;
//        std::copy(
//                make_choose_polynomial_iterator(F.begin()),
//                make_choose_polynomial_iterator(F.end()),
//                std::back_inserter(result));
        using namespace std::tr1::placeholders;
        std::transform(F.begin(), F.end(), std::back_inserter(result),
                std::tr1::bind(&PointPolyPair::second, _1));
        return result;
    }

    typedef typename PointPolyMap::value_type PointPolyPair;

    typedef std::tr1::function<Point<Dim, OrderPolicy> const & (PointPolyPair const &)>
        ChoosePointFunctor;

    template<typename It>
    boost::transform_iterator<ChoosePointFunctor, It>
    make_choose_point_iterator(It it) {
        using std::tr1::bind;
        using namespace std::tr1::placeholders;
        ChoosePointFunctor choosePoint(std::tr1::bind(&PointPolyPair::first, _1));
        return boost::make_transform_iterator(it, choosePoint);
    }
};

//template<typename T, template <typename PointImpl> class OrderPolicy>
//const typename BMSAlgorithm<T, OrderPolicy>::CoefT
//BMSAlgorithm<T, OrderPolicy>::ZERO =
//    CoefficientTraits<typename BMSAlgorithm<T, OrderPolicy>::CoefT>::addId();

template<typename SeqT,typename PolynomialT, template <typename PointImpl> class OrderPolicy>
BMSAlgorithm<SeqT, PolynomialT, OrderPolicy>
make_bmsalgorithm(
        SeqT const & seq,
        Point<BMSAlgorithm<SeqT, PolynomialT, OrderPolicy>::Dim, OrderPolicy>
            const & seqLen) {
    return BMSAlgorithm<SeqT, PolynomialT, OrderPolicy>(seq, seqLen);
}

template<typename SeqT,typename PolynomialT>
BMSAlgorithm<SeqT, PolynomialT>
make_bmsalgorithm(
        SeqT const & seq,
        Point<BMSAlgorithm<SeqT, PolynomialT>::Dim>
            const & seqLen) {
    return BMSAlgorithm<SeqT, PolynomialT>(seq, seqLen);
}

template<typename SeqT>
BMSAlgorithm<SeqT>
make_bmsalgorithm(
        SeqT const & seq,
        Point<BMSAlgorithm<SeqT>::Dim>
            const & seqLen) {
    return BMSAlgorithm<SeqT>(seq, seqLen);
}

bool defPointsOrder(Point<2> const & lhs, Point<2> const & rhs) {
    return lhs[0] < rhs[0] && lhs[1] > rhs[1];
}

template<typename CoefT>
std::list< Point<2> > getOrderedDefPoints(
        typename BMSAlgorithm<
            typename MVPolyType<2, CoefT>::type
            >::PointPolyMap const & ppmap) {
    typedef typename BMSAlgorithm<
                typename MVPolyType<2, CoefT>::type
            >::PointPolyMap::value_type map_value_type;
    using namespace std::tr1::placeholders;
    std::list< Point<2> > result;
    transform(
            ppmap.begin(),
            ppmap.end(),
            std::back_inserter(result),
            std::tr1::bind( &map_value_type::first, _1 )
    );
    result.sort(defPointsOrder);
    return result;
}

void printOrderedDefPoints(std::list< Point<2> > const & points) {
    if (points.empty())
        return;
    using std::cout;
    using std::endl;
    using std::setw;

    cout << setw(4) << "";
    for (int i = 0; i != (*points.begin())[1] + 1; ++i) {
        cout << setw(2) << i;
    }
    cout << endl;
    cout << setw(4) << "";
    for (int i = 0; i != (*points.begin())[1] + 1; ++i) {
        cout << "__";
    }
    cout << endl;

    std::list< Point<2> >::const_iterator it = points.begin();
    const int end = (*--points.end())[0] + 1;
    for (int i = 0; i != end; ++i) {
        cout << setw(2) << i << " | ";
        if ((*it)[0] == i)
            cout << std::string((*it++)[1] * 2 + 1, '~') << "+";
        cout << endl;
    }
}

} // namespace mv_poly

#endif /* BMSA_HPP_ */
