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

#include <list>
#include <map>

#include <tr1/functional>

#include <boost/iterator/transform_iterator.hpp>

#include "mv_poly.hpp"
#include "Point.hpp"
#include "CoefficientTraits.hpp"

template<typename PolynomialT>
class BMSAlgorithm {
    static const int Dim = PolynomialT::VAR_CNT;

    typedef std::list< Point<Dim> > PointCollection;

    typedef std::list< PolynomialT > PolynomialCollection;

    typedef std::map< Point<Dim>, PolynomialT > PointPolyMap;

    PointPolyMap F, G;

    typedef typename PolynomialT::CoefT CoefT;

    static const CoefT ZERO = CoefT();

    typedef std::map< Point<Dim>, CoefT > PointCoefMap;

    PolynomialT const & seq;

    Point<Dim> seqLen;

public:

    BMSAlgorithm(PolynomialT const & seq_, Point<Dim> const & seqLen_) :
        seq(seq_), seqLen(seqLen_)  {}

    PolynomialCollection computeMinimalSet() {
        using std::tr1::bind;
        using namespace std::tr1::placeholders;
        PointCollection deltaPoints, sigmaPoints;
        PointCoefMap discrepancies;
        typedef typename PointPolyMap::value_type PointPolyPair;
        std::tr1::function<
            Point<Dim> const & (PointPolyPair const &)
            > choosePoint(std::tr1::bind(&PointPolyPair::first, _1));
        for (Point<Dim> k; k < seqLen; ++k) {
            PointPolyMap NewF, NewG;
            for (typename PointPolyMap::const_iterator fIt = F.begin();
                    fIt != F.end(); ++fIt) {
                Point<Dim> const & degF = fIt->first;
                PolynomialT const & f = fIt->second;
                if (byCoordinateLess(degF, k)) {
                    CoefT b = conv(f, seq, degF, k);
                    Point<Dim> d = k - degF;
                    deltaPoints.push_back(d);
                    if (b != ZERO &&
                            !byCoordinateLessThenAny(
                                    d,
                                    boost::make_transform_iterator(
                                            F.begin(), choosePoint),
                                    boost::make_transform_iterator(
                                            F.end(), choosePoint))) {

                    }
                }
            }
            deltaPoints = getPartialMaximums(deltaPoints);
            sigmaPoints = getConjugatePointCollection(deltaPoints);
            for (typename PointCollection::const_iterator
                    cIt = deltaPoints.begin();
                    cIt != deltaPoints.end(); ++cIt) {
                typename PointPolyMap::const_iterator tmpIt = G.find(*cIt);
                if (tmpIt != G.end())
                    NewG.insert(*tmpIt);
                else {
                    Point<Dim> s = k - *cIt;
                    newG[*cIt] =
                            CoefficientTraits<CoefT>::inverse(discrepancies[s])
                            * F[s];
                }
            }
            for (typename PointCollection::const_iterator
                    sIt = sigmaPoints.begin();
                    sIt != sigmaPoints.end(); ++sIt) {

            }
        }

        return getPolynomialList();

    }

private:
    PolynomialCollection getPolynomialList() {
        PolynomialCollection result;
        // TODO: method stub
        return result;
    }
};

#endif /* BMSA_HPP_ */
