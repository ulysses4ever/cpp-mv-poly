/*
 * BasicExamples.hpp
 *
 *  Created on: 21.06.2011
 *      Author: ulysses
 */

#ifndef BASICEXAMPLES_HPP_
#define BASICEXAMPLES_HPP_

#include <algorithm>
#include <deque>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <sstream>
#include <set>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include <cmath>

#include <tr1/functional>
#include <boost/foreach.hpp>

#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>

#include "mv_poly.hpp"
#include "Point.hpp"
#include "bmsa.hpp"
#include "NtlUtilities.hpp"
#include "NtlPolynomials.hpp"

using namespace mv_poly;

using namespace std;

typedef MVPolyType<2, NTL::GF2E>::ResultT PolyT;

PolyT sakataEtAl05Seq(NTL::GF2E const & x) {
    // prepare raw coefficient arrays for reference sequence
    using NTL::GF2E;
    using NTL::power;
    GF2E coefs0[] = {
            power(x, 9), GF2E(), power(x, 9), power(x, 6), power(x, 5),
            power(x, 6)
    };
    GF2E coefs1[] = {
            power(x, 14), power(x, 9), power(x, 11), power(x, 4), power(x, 7)
    };
    GF2E coefs2[] = {
            power(x, 5), power(x, 14), GF2E(), power(x, 7)
    };
    GF2E coefs3[] = {
            power(x, 7), power(x, 12), power(x, 12)
    };
    GF2E coefs4[] = {
            power(x, 2), power(x, 5)
    };
    GF2E coefs5[] = {
            power(x, 5), power(x, 5)
    };
    GF2E coefs6[] = {
            power(x, 0)
    };
    GF2E const * coefs[] = {coefs0, coefs1, coefs2, coefs3, coefs4, coefs5, coefs6};
    size_t coefsLens[] = {
            ARR_LEN(coefs0), ARR_LEN(coefs1), ARR_LEN(coefs2),
            ARR_LEN(coefs3), ARR_LEN(coefs4), ARR_LEN(coefs5), ARR_LEN(coefs6)
    };

    // load reference sequence
    return load_coefs<PolyT>(coefs, ARR_LEN(coefs), coefsLens);
}

void testField() {
    //init field
    typedef NTL::GF2 PrimeField;
    typedef typename NTLPrimeFieldTtraits<PrimeField>::ExtField ExtField;
    size_t field_size = initExtendedField<PrimeField>("[1 1 0 0 1]");
    //cout << field_size;

    // primitive element
    ExtField x = getPrimitive<ExtField>();

    printFieldInPowers(x, field_size);
    MVPolyType<2, ExtField>::type u = sakataEtAl05Seq(x);
    Point<2> pt;
    pt[0] = 1; pt[1] = 4;
    cout << u[pt] << endl;
    ExtField a = u[pt];
    pt[1] = 1;
    cout << u[pt] << endl;
    cout << a + u[pt] << endl;
}

void testCoxEtAl05() {

    //init field
    typedef NTL::GF2 PrimeField;
    typedef typename NTLPrimeFieldTtraits<PrimeField>::ExtField ExtField;
    // size_t field_size =
    initExtendedField<PrimeField>("[1 1 1]");
    //cout << field_size;

    // primitive element
    ExtField x = getPrimitive<ExtField>();

    // prepare coefs for reference sequence
    ExtField coefs0[] = {
            ExtField::zero(), x
    };
    ExtField coefs1[] = {
            power(x, 2), x
    };
    ExtField coefs2[] = {
            x
    };
    ExtField const * coefs[] = {coefs0, coefs1, coefs2};
    size_t coefsLens[] = {
            ARR_LEN(coefs0), ARR_LEN(coefs1), ARR_LEN(coefs2)
    };

    // load reference sequence
    const int dim = 2;
    typedef MVPolyType<dim, ExtField>::ResultT PolyT;
    PolyT u = load_coefs<PolyT>(coefs, ARR_LEN(coefs), coefsLens);

    // need to manually set sequence length
    Point<2> ulen;
    ulen[0] = 0; ulen[1] = 2;

    // print reference sequence
    //cout << u << endl;

    // run algorithm
    BMSAlgorithm< PolyT > alg(u, ulen);
    BMSAlgorithm< PolyT >::PolynomialCollection minset = alg.computeMinimalSet();
    BMSAlgorithm< PolyT >::PointPolyMap F = alg.getF();

    // print algo result
//    cout << endl << "copmuted res:" << endl;
//    copy(minset.begin(), minset.end(),
//                std::ostream_iterator<PolyT>(cout, "\n"));
    // [[[1] [1]] [[1 1]]]  ~ x_2 + (a^2)x_1 + 1
    // [[[]] [[1 1]] [[1]]] ~ (x_1)^2 + (a^2)x_1

    BOOST_FOREACH(PolyT const & f, minset) {
        cout << makePowerPrinter< GradedAntilexMonomialOrder >(f, x) << endl;
    }


}

void testSakataEtAl95() {
    std::istringstream is("[1 1 0 0 1]"); // deg = 4
    //const size_t field_size = 16;
    NTL::GF2X irr;
    is >> irr;
//    cout << irr << endl;
    NTL::GF2E::init(irr);

    is.str("[0 1]");
    NTL::GF2E x;
    using std::cout;
    using std::endl;
//    std::cout << "No-arg ctor of GF2: " << x << endl;
    is >> x;
//    cout << x << endl;

    // load reference sequence
    typedef MVPolyType<2, NTL::GF2E>::ResultT PolyT;
    PolyT u = sakataEtAl05Seq(x); //load_coefs<PolyT>(coefs, ARR_LEN(coefs), coefsLens);
    Point<2> ulen;
    //ulen[0] = 4; ulen[1] = 2;
    ulen[0] = 5; ulen[1] = 1;
    //u.setCoefs(stor);
    cout << u << endl;

    PolyT expected_res3;
    std::list<PolyT> expected_res;
    {   // expected res

        PolyT expected;
        std::ostringstream oss;
        std::istringstream iss;

        oss.str("");
        oss
                << '['
                    << '['
                        << power(x, 11) << ' ' << power(x, 13) << ' ' << power(x, 0)
                    << ']'
                    << '['
                        << power(x, 13) << ' ' << power(x, 10)
                    << ']'
                << ']';
        iss.str(oss.str());
        iss >> expected;
        expected_res.push_back(expected);

        expected_res3 = expected;

        oss.str("");
        oss
                << '['
                    << '['
                        << power(x, 4) << ' ' << power(x, 1)
                    << ']'
                    << '['
                        << power(x, 7) << ' ' << power(x, 4)
                    << ']'
                    << '['
                        << power(x, 3) << ' ' << power(x, 0)
                    << ']'
                << ']';
        iss.str(oss.str());
        iss >> expected;
        expected_res.push_back(expected);

        oss.str("");
        oss
                << '['
                    << '['
                        << power(x, 1) << ' ' << power(x, 11)
                    << ']'
                    << '['
                        << power(x, 1) << ' ' << power(x, 0)
                    << ']'
                    << '['
                        << power(x, 4) << ' ' << power(x, 9)
                    << ']'
                    << '['
                        << power(x, 3)
                    << ']'
                    << '['
                        << power(x, 0)
                    << ']'
                << ']';
        iss.str(oss.str());
        iss >> expected;
        expected_res.push_back(expected);
    }

//    cout <<  endl <<  endl;
    //    cout << CoefficientTraits<GF2E>::addId() << endl;
//    cout << CoefficientTraits<GF2E>::multId() << endl;
//    cout << CoefficientTraits<NTL::GF2>::addId() << endl;
//    cout << CoefficientTraits<NTL::GF2>::multId() << endl;
//    GF2E zero;
//    NTL::clear(zero);
//    cout << zero << endl;

    BMSAlgorithm< PolyT > alg(u, ulen);
    BMSAlgorithm< PolyT >::PolynomialCollection minset = alg.computeMinimalSet();
    BMSAlgorithm< PolyT >::PointPolyMap F = alg.getF();

    //assert(expected_res == minset);
    cout << "expected res:" << endl;
//    copy(expected_res.begin(), expected_res.end(),
//            std::ostream_iterator<PolyT>(cout, "\n"));
//    std::vector<int> expectes_res_powers;
//    transform(expected_res.begin(), expected_res.end(),
//            back_inserter(expectes_res_powers),
//            unmap_each(powers, field_size));
//    copy(expected_res.begin(), expected_res.end(),
//            std::ostream_iterator<PolyT>(cout, "\n"));

    cout << endl << "copmuted res:" << endl;
//    copy(minset.begin(), minset.end(),
//                std::ostream_iterator<PolyT>(cout, "\n"));
    for (std::list<PolyT>::const_iterator it = minset.begin();
            it != minset.end(); ++it) {
        cout << makePowerPrinter< GradedAntilexMonomialOrder >(*it, x)
                << endl;
    }

    // print all convolutions u * expected - should be 0's
//    {
//        PolyT poly =  minset.back();
//        Point<2> deg_poly;
//        deg_poly[0] = 4; deg_poly[1] = 0;
//        for (Point<2> k = deg_poly; k < ulen; ++k) {
//            if (byCoordinateLess(deg_poly, k))
//                cout << "k = " << k << ", conv(f, u, k) = "
//                    << conv(expected_res3, u, deg_poly, k) << endl;
//        }
//        cout << endl;
//    }

    BMSAlgorithm< PolyT >::PolynomialCollection::const_iterator it = minset.begin();
    list<PolyT>::const_iterator it_exp = expected_res.begin();
    assert(*it++ == *it_exp++);
    assert(*it++ == *it_exp++);
//    assert(*it++ == *it_exp++);
    //assert(minset2 == expected_res);
//    copy(F.begin(), F.end(), std::ostream_iterator<
//            BMSAlgorithm< PolyT >::PointPolyMap::value_type>(cout, "\n"));
//    cout << *it << endl;
//    cout << *it_exp << endl;
    //printOrderedDefPoints(getOrderedDefPoints<NTL::GF2E>(F));


//    copy(minset2.begin(), minset2.end(), std::ostream_iterator<PolyT>(cout, "\n"));
}

void miscTesting() {
//    typedef MVPolyType<2, NTL::GF2>::ResultT PolyT;
//    PolyT u("[[0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]");
//    Point<2> pt;
//    pt[0] = 4; pt[1] = 1;
//    BMSAlgorithm< PolyT > alg(u, pt);
//    BMSAlgorithm< PolyT >::PolynomialCollection minset2 = alg.computeMinimalSet();
//    BMSAlgorithm< PolyT >::PointPolyMap F = alg.getF();
//    copy(minset2.begin(), minset2.end(), std::ostream_iterator<PolyT>(cout, "\n"));
//
//    Point<3> ptt;
//    ptt[0] = 5; ptt[1] = 0; ptt[2] = 1;
//    typedef MVPolyType<3, NTL::GF2>::ResultT PolyT3;
//    PolyT3 v(
//            "[[[1 1 1 1 0 0] [0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]"//x = 0
//            "[[1 1 0 1 1] [1 0 1 1] [0 1 1] [1 1] [1] [0]]" // x = 1
//            "[[0 1 0 0] [0 0 1] [0 0] [1] [0]]" // x = 2
//            "[[1 1 0] [1 0] [0] [1]] [[1 1] [0] [1]] [[1] [1]] [[0]]]" // x = 3-6
//            );
//    BMSAlgorithm< PolyT3 > alg3(v, ptt);
//    BMSAlgorithm< PolyT3 >::PolynomialCollection minset = alg3.computeMinimalSet();
//    //BMSAlgorithm< PolyT3 >::PointPolyMap F3 = alg3.getF();
//
//    copy(minset.begin(), minset.end(), std::ostream_iterator<
//            BMSAlgorithm< PolyT3 >::PolynomialCollection::value_type>(cout, "\n"));

    //    Point<3> degf;
//    degf[0] = 0; degf[1] = 0; degf[2] = 3;
    //PolyT3 const & f = F3[degf];
//    PolyT3 f("[[[1] [1 1]]]");

//    for (Point<3> i; i < ptt; ++i) {
//        if (byCoordinateLess(degf, i))
//            cout << "i: " << i << " discr = " << conv(f, v, degf, i) << endl;
//    }
//    copy(F3.begin(), F3.end(), std::ostream_iterator<
//            BMSAlgorithm< PolyT3 >::PointPolyMap::value_type>(cout, "\n"));
//
//    PolyT3 f("[[[1] [1]] [[1]]]");
//    ptt[0] = 0; ptt[1] = 1; ptt[2] = 1;
//    degf[0] = 0; degf[1] = 1; degf[2] = 0;
//    cout << conv(f, v, degf, ptt) << endl;

//    const int Dim = 3;
//    Point<Dim> upperPoint;
//    upperPoint[0] = 1 + 2;
//    std::list< Point<Dim> > approxSigmaSet;
//    for (Point<Dim> i; i < upperPoint; ++i) {
////        if (! byCoordinateLessThenAny(i, points))
////            approxSigmaSet.push_back(i);
//        //cout << "i: " << i << endl;
//    }
//
//    Point<4> pt, up;
//    pt[0] = 4;
//    up[0] = 5;
//    while(pt++ < up) {
//        cout << pt << endl;
//    }
//    ptt[0] = 2; ptt[1] = 0; ptt[2] = 1;
//    PolyT3 f("[[[1 0 1] [0 1]] [[1 0] [1 1]]]");
//    f <<= ptt;
//    (( 0 1 ), [[1 1] [1]])
//    (( 2 0 ), [[0] [0] [1]])

//    k = ( 2 1 )
//    failed: [[1 1] [1]]
//    span: ( 2 0 )

//    typedef std::map<Point<2>, PolyT> PointPolyMap;
//    PointPolyMap m;
//    pt[0] = 0; pt[1] = 1;
//    m.insert(std::make_pair(pt, PolyT("[[0 1]]")));
//    pt[0] = 2; pt[1] = 0;
//    m.insert(std::make_pair(pt, PolyT("[[0] [0] [1]]")));
//
//    pt[0] = 2; pt[1] = 0;
//    cout << m[pt] << endl;
//    copy(m.begin(), m.end(),
//            std::ostream_iterator<PointPolyMap::value_type>(cout, "\n"));

//    Point<2> pt4, pt5;
//    pt4[0] = 0; pt4[1] = 1;
//    pt5[0] = 2; pt5[1] = 0;
//    cout << totalLess(pt4, pt5) << endl;
//    cout << totalLess(pt5, pt4) << endl;
//    (( 2 0 ), [[0] [0] [1]])
//    (( 0 1 ), [[0 1]])

//    pt[0] = 0; pt[1] = 0;
//    PolyT f("[[0] [0] [1]]");
//    PolyT g = f << pt;
//    cout << g << endl;

//    MVPolyType<2, int>::ResultT p("[[1 0 1] [1 1]]"), q("[[2 3] [0 2] [3]]");
//    MVPolyType<2, int>::ResultT pPlusQ = p + q;// "[[3 3 1]] [1 3] [3]";
//    cout << pPlusQ << endl;

//    MVPolyType<2, int>::ResultT p("[[1 0 1] [1 1]]"), q(p);
//    q += p;
//    p *= 2;
//    cout << p << endl;
//    cout << q << endl;
//    cout << (p == q) << endl;
//    cout << p << endl;
//    cout << q << endl;

//    MVPolyType<2, int>::ResultT p("[[1 0 1] [1 1]]");
//    Point<2> pt;
//    pt[0] = 0; pt[1] = 1;
//    cout << (p << pt) << endl;
//    pt[0] = 1; pt[1] = 0;
//    cout << (p << pt) << endl;

//    Point<2> pt;
//    std::list<Point<2> > s, sn, sig;
//    pt[0] = 0; pt[1] = 1;
//    s.push_back(pt);
//    pt[0] = 2; pt[1] = 0;
//    s.push_back(pt);
//    pt[0] = 1; pt[1] = 0;
//    s.push_back(pt);
//    sn = getPartialMaximums(s);
//
//    sig = getConjugatePointCollection(sn);
    //copy(sn.begin(), sn.end(), ostream_iterator<Point<2> >(cout, "\n"));

    //std::copy(sn.begin(), sn.end(), ostream_iterator<Point<2> >(cout, "\n"));
    //std::cout << std::endl;
//    std::copy(sig.begin(), sig.end(), ostream_iterator<Point<2> >(cout, "\n"));

//    MVPolyType<2, NTL::GF2>::ResultT f, u;
//    loadPolyFromString(u, "[[0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]");
//    loadPolyFromString(f, "[[1 1] [1]]");
//    Point<2> degf;
//    degf[0] = 0;
//    degf[1] = 1;
//    Point<2> m;
//    m[0] = 0;
//    m[1] = 2;
//    cout << conv(f, u, degf, m) << endl;
//    m[0] = 2;
//    m[1] = 1;
//    cout << conv(f, u, degf, m) << endl;
//    for ( ; totalLess(i, pt); increase(i)) {
//        cout << p[i] << " ";
//    }
//    typedef MVPolyType<1, int>::ResultT Poly1;
//    Poly1 p1;
//    loadPolyFromString(p1, "[1 2 3]");
    //cout << p1 << endl;

//    MVPolyType<2, int>::ResultT p2;
//    loadPolyFromString(p2, "[[1 2 3] [3 2 1] [1]]");
//    cout  << endl << p2 << endl;

//    MVPolyType<3, int>::ResultT p3;
//    loadPolyFromString(p3, "[[[1 2] [3]] [[3] [2 42]] [[1]]]");
//    os.str("");
//    os << p2;
//    s = "[[[1 2] [3]] [[3] [2 1]] [[1]]]";
//    Point<3> pt;
//    pt[0] = 1; pt[1] = 1; pt[2] = 1;
//    cout << p3[pt];

//    Point<3> pt1, pt2, pt3;
//    pt1[0] = 3; pt1[1] = 1; pt1[2] = 2;
//    pt2[0] = 2; pt2[1] = 1; pt2[2] = 0;
//    pt3[0] = 2; pt3[1] = 2; pt3[2] = 0;
//    cout << pt1 << endl;
//    cout << byCoordinateLess(pt2, pt1) << endl;  // pt1 <_p pt2
//    cout << !byCoordinateLess(pt1, pt2) << endl; // pt2 \not <_p pt1 as pt1 <_p pt2 (see above)
//    cout << !byCoordinateLess(pt1, pt3) << endl; // pt1 and pt3 incomparable

//    Point<2> pt4, pt5;
//    pt4[0] = 1; pt4[1] = 0;
//    pt5[0] = 0; pt4[1] = 1;
//    cout << totalLess(pt4, pt5) << endl;

//    Point<3> pt;
//    pt[0] = 0; pt[1] = 0; pt[2] = 0;
//    for (int i = 0; i < 12; ++i) {
//        cout << pt << endl;
//        increase(pt);
//    }
}

/// ---------------------- enable_if exercises ---------------------------

//template<typename T>
//struct A {
//
//    template<typename S>
//    typename boost::enable_if<boost::is_arithmetic <S>, T>::type f(S s) {}
//
//    template<typename S>
//    typename boost::enable_if<boost::is_array<S>, T>::type f(S s) {}
//};

//template<typename S>
//typename boost::enable_if_c<true, int>::type f(S s) {}
//
//template<typename S>
//typename boost::disable_if_c<false, int>::type f(S s) {}

//
//template <typename T>
//void pass(T t, typename boost::enable_if< boost::is_arithmetic <T> >::type * = 0) //boost::is_arithmetic <T>
//{
//// . . .
//}
//template <typename T>
//void pass(T t, typename boost::disable_if< boost::is_integral<T> >::type * = 0)
//{
//// ...
//}

//template <typename T>
//void pass_1(T t, typename boost::enable_if< boost::mpl::true_ >::type * = 0) //boost::is_arithmetic <T>
//{
//// . . .
//}
//template <typename T>
//void pass_1(T t, typename boost::disable_if< boost::mpl::false_ >::type * = 0)
//{
//// ...
//}


#endif /* BASICEXAMPLES_HPP_ */
