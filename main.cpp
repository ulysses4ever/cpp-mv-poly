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
#include <sstream>
#include <string>
#include <typeinfo>

//#include <boost/test/>

#include <NTL/GF2.h>

#include "mv_poly.hpp"
#include "Point.hpp"

using std::cout;
using std::endl;
using std::string;

template<typename T>
class Bar {
public:
    //static const int =
};

int main() {
    MVPolyType<2, NTL::GF2>::ResultT p;
    loadPolyFromString(p, "[[0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]");
    Point<2> i, pt;
    pt[0] = 4; pt[1] = 1;
    for ( ; totalLess(i, pt); increase(i)) {
        cout << p[i] << " ";
    }
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
