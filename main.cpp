/*
 * main.cpp
 *
 *  Created on: 06.12.2009
 *      Author: ulysses
 */

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

#include <tr1/functional>

#include <NTL/GF2.h>

#include "mv_poly.hpp"
#include "Point.hpp"
#include "bmsa.hpp"

using std::cout;
using std::endl;
using std::string;
using std::copy;
using std::ostream_iterator;

int main() {
//    typedef MVPolyType<2, NTL::GF2>::ResultT PolyT;
//    PolyT u("[[0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]");
//    Point<2> pt;
//    pt[0] = 4; pt[1] = 1;
//    BMSAlgorithm< PolyT > alg(u, pt);
//    BMSAlgorithm< PolyT >::PolynomialCollection minset = alg.computeMinimalSet();
//    BMSAlgorithm< PolyT >::PointPolyMap F = alg.getF();
//    copy(F.begin(), F.end(), std::ostream_iterator<
//            BMSAlgorithm< PolyT >::PointPolyMap::value_type>(cout, "\n"));
    //copy(minset.begin(), minset.end(), std::ostream_iterator<PolyT>(cout, "\n"));

    Point<3> ptt;
    ptt[0] = 5; ptt[1] = 0; ptt[2] = 1;
    typedef MVPolyType<3, NTL::GF2>::ResultT PolyT3;
    PolyT3 v(
            "[[[1 1 1 1 0 0] [0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]"//x = 0
            "[[1 1 0 1 1] [1 0 1 1] [0 1 1] [1 1] [1] [0]]" // x = 1
            "[[0 1 0 0] [0 0 1] [0 0] [1] [0]]" // x = 2
            "[[1 1 0] [1 0] [0] [1]] [[1 1] [0] [1]] [[1] [1]] [[0]]]" // x = 3-6
            );
    BMSAlgorithm< PolyT3 > alg3(v, ptt);
    BMSAlgorithm< PolyT3 >::PolynomialCollection minset = alg3.computeMinimalSet();
    BMSAlgorithm< PolyT3 >::PointPolyMap F3 = alg3.getF();
    copy(F3.begin(), F3.end(), std::ostream_iterator<
            BMSAlgorithm< PolyT3 >::PointPolyMap::value_type>(cout, "\n"));

    ptt[0] = 2; ptt[1] = 0; ptt[2] = 1;
    PolyT3 f("[[[1 0 1] [0 1]] [[1 0] [1 1]]]");
    f <<= ptt;
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
