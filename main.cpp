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

#include <boost/regex.hpp>
//#include <boost/test/>

#include "mv_poly.hpp"
#include "Point.hpp"

using std::cout;
using std::endl;
using std::string;

int main() {
//    typedef MVPolyType<1, int>::ResultT Poly1;
//    Poly1 p1;
//    loadPolyFromString(p1, "[1 2 3]");
//    cout << p1 << endl;

//    MVPolyType<2, int>::ResultT p2;
//    loadPolyFromString(p2, "[[1 2 3] [3 2 1] [1]]");
//    cout  << endl << p2 << endl;

//    MVPolyType<3, int>::ResultT p3;
//    loadPolyFromString(p2, "[[[1 2] [3]] [[3] [2 1]] [[1]]]");
//    os.str("");
//    os << p2;
//    s = "[[[1 2] [3]] [[3] [2 1]] [[1]]]";
    //Point<5> p;

//    Point<3> pt1, pt2, pt3;
//    pt1[0] = 3; pt1[1] = 1; pt1[2] = 2;
//    pt2[0] = 2; pt2[1] = 1; pt2[2] = 0;
//    pt3[0] = 2; pt3[1] = 2; pt3[2] = 0;
//    cout << pt1 << endl;
//    cout << byCoordinateLess(pt2, pt1) << endl;  // pt1 <_p pt2
//    cout << !byCoordinateLess(pt1, pt2) << endl; // pt2 \not <_p pt1 as pt1 <_p pt2 (see above)
//    cout << !byCoordinateLess(pt1, pt3) << endl; // pt1 and pt3 incomparable

    Point<2> pt4, pt5;
    pt4[0] = 1; pt4[1] = 0;
    pt5[0] = 0; pt4[1] = 1;
    cout << totalLess(pt4, pt5) << endl;

}
