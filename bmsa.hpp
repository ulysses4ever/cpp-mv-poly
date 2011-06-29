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

#include <tr1/functional>

#include <boost/iterator/transform_iterator.hpp>

#include "mv_poly.hpp"
#include "Utilities.hpp"
#include "Point.hpp"
#include "CoefficientTraits.hpp"

namespace mv_poly {

template<
    typename PolynomialT,
    template <typename PointImpl> class OrderPolicy = GradedAntilexMonomialOrder
>
class BMSAlgorithm {
public:
    static const int Dim = PolynomialT::VAR_CNT;

    typedef std::map< Point<Dim, OrderPolicy>, PolynomialT > PointPolyMap;

private:
    typedef std::list< Point<Dim, OrderPolicy> > PointCollection;

    typedef typename PolynomialT::CoefT CoefT;

    typedef std::map< Point<Dim, OrderPolicy>, CoefT > PointCoefMap;

    PointPolyMap G;

    const CoefT ZERO;

    Point<Dim, OrderPolicy> seqLen;

    PointCollection oldDeltaPoints;

protected:

    PointPolyMap F;

    PolynomialT const & seq;

    void infoUpdate(Point<Dim> const & k) {
        using std::tr1::bind;
        using std::tr1::cref;
        using namespace std::tr1::placeholders;
        using std::cout;
        using std::endl;

//            cout << "k = " << k << endl;
            PointPolyMap newF, newG;
            PointCollection deltaPoints, sigmaPoints;
            PointCoefMap discr; // discrepancies
            //searching for candidates to form new deltaPoints
            for (typename PointPolyMap::const_iterator fIt = F.begin();
                    fIt != F.end(); ++fIt) {
                Point<Dim> const & degF = fIt->first;
                PolynomialT const & f = fIt->second;
                if (byCoordinateLess(degF, k)) {
                    CoefT b = conv(f, seq, degF, k);
                    discr[degF] = b;
                    Point<Dim> c = k - degF;
                    if (b != ZERO) {
//                        cout << "failed: " << f << endl;
//                        cout << "span: " << c << endl;
                    }
                    if (b != ZERO &&
                            !byCoordinateLessThenAny(
                                    c,
                                    make_choose_point_iterator(G.begin()),
                                    make_choose_point_iterator(G.end()))) {
                        deltaPoints.push_back(c);
//                        cout << "fallen: " << f << endl;
                    }
                }
            }
            deltaPoints.splice(deltaPoints.end(), oldDeltaPoints);
            deltaPoints = getPartialMaximums(deltaPoints);
            sigmaPoints = getConjugatePointCollection(deltaPoints);

//            cout << "Delta-points:" << endl;
//            copy(deltaPoints.begin(), deltaPoints.end(),
//                    std::ostream_iterator< Point<Dim> >(cout, "\n"));
//            cout << "Sigma-points:" << endl;
//            copy(sigmaPoints.begin(), sigmaPoints.end(),
//                    std::ostream_iterator< Point<Dim> >(cout, "\n"));

            // forming new G
            for (typename PointCollection::const_iterator
                    cIt = deltaPoints.begin();
                    cIt != deltaPoints.end(); ++cIt) {
                typename PointPolyMap::const_iterator tmpIt = G.find(*cIt);
                if (tmpIt != G.end())
                    newG.insert(*tmpIt);
                else {
                    Point<Dim> s = k - *cIt;
                    newG[*cIt] =
                            CoefficientTraits<CoefT>::multInverse(discr[s])
                            * F[s];
                }
            }

            // forming new F
            for (typename PointCollection::const_iterator
                    tIt = sigmaPoints.begin();
                    tIt != sigmaPoints.end(); ++tIt) {
                Point<Dim> const & t = *tIt;
                // next find_if will find something for shure — I GUARANTEE IT
                Point<Dim> s = *(std::find_if(
                        make_choose_point_iterator(F.begin()),
                        make_choose_point_iterator(F.end()),
                        bind(&byCoordinateLess<Dim, OrderPolicy>, _1, cref(t))));
                Point<Dim> u = t - s; // valid Point as true == byCoordinateLess(s, t)
//                cout << "t: " << t << " ";
//                cout << "s: " << s << " ";
//                cout << "u: " << u;
//                cout << endl;
                //cout << "F[s] << u: " << (F[s] << u) << endl;
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
                        Point<Dim> const & c = cIt->first;
                        newF[t] = (F[s] << u) - // Berlekamp formula
                            (discr[s] * G[c] << (c - (k - t)));
                    }
                }
                if (!notJustIncreaseDegree) {
                    newF[t] = F[s] << u;
//                    cout << "Just increased degree on " << u << endl;
                }
            }
            F = newF;
            G = newG;
            oldDeltaPoints.clear();
            oldDeltaPoints.splice(oldDeltaPoints.end(), deltaPoints);

//            cout << "F: " << endl;
//            copy(F.begin(), F.end(),
//                    std::ostream_iterator<typename PointPolyMap::value_type>(cout, "\n"));
//            cout << "G: " << endl;
//            copy(G.begin(), G.end(),
//                    std::ostream_iterator<typename PointPolyMap::value_type>(cout, "\n"));
//            cout << endl;
    }

public:

    typedef std::list< PolynomialT > PolynomialCollection;

    BMSAlgorithm(PolynomialT const & seq_,
            Point<Dim, OrderPolicy> const & seqLen_) :
                ZERO(CoefficientTraits<CoefT>::addId()),
                seqLen(seqLen_), seq(seq_) {
        F.insert(std::make_pair(Point<Dim>(), PolynomialT::getId()));
    }

    PolynomialCollection computeMinimalSet() {
        oldDeltaPoints.clear();
        // scanning input sequense step-by-step, following monomial order
        for (Point<Dim> k; k < seqLen; ++k) {
            infoUpdate(k);
        }
        return getPolynomialList();

    }

    PointPolyMap const & getF() {
        return F;
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

    typedef std::tr1::function<Point<Dim> const & (PointPolyPair const &)>
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

template<typename PolynomialT, template <typename PointImpl> class OrderPolicy>
BMSAlgorithm<PolynomialT, OrderPolicy>
make_bmsalgorithm(PolynomialT const & seq,
        Point<BMSAlgorithm<PolynomialT, OrderPolicy>::Dim, OrderPolicy> const & seqLen) {
    return BMSAlgorithm<PolynomialT, OrderPolicy>(seq, seqLen);
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
