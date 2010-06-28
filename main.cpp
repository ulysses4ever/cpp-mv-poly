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

using std::cout;
using std::endl;
using std::string;

int main() {
    typedef MVPolyType<1, int>::ResultT Poly1;
    Poly1 p1;
    //std::string s("[[1 2 3][4 5 6]]");
    //std::istringstream iss();
    //std::istringstream iss("[4 5 6]");
    string s("[1 2 3]");
    loadPolyFromString(p1, s);
    std::cout << "Число переменных Poly1: " << Poly1::VAR_CNT << std::endl;
    std::cout << "Тип Poly1: " << typeid(p1).name() << std::endl;
    std::cout << "p1: " << p1 << std::endl;

}
