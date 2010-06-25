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

template<typename T>
class VarCnt {
public:
    static const int result = 0;
};

template<typename T>
class Polynomial {
public:
    typedef T                   CoefT;
    typedef std::deque<CoefT>   StorageT;
    static const int VAR_CNT = 1 + VarCnt<T>::result;

    template<typename S>
    friend
    void loadPolyFromString(Polynomial<S> & p, std::string const & s );

    Polynomial() {}

    Polynomial(std::string const & s) {
        loadPolyFromString(*this, s );
    }

    void setCoefs(StorageT const & data)  { this->data = data; }
    StorageT const & getCoefs() const     { return data; }

private:
    StorageT data;
};

template<typename T>
class VarCnt< Polynomial<T> > {
public:
    static const int result = Polynomial<T>::VAR_CNT;
};

template<typename T>
void loadPolyFromString( Polynomial<T> & p, std::string const & s ) {
    std::istringstream is(s);
    if (is.peek() == '[')
        is.get();
    else {
        is.setstate(std::ios::failbit);
        return;
    }
    typename Polynomial<T>::StorageT tempStorage;

    typename Polynomial<T>::CoefT cf;
    while ( is >> cf ) {
        tempStorage.push_back(cf);
    }

    p.setCoefs( tempStorage );

}

/*
using std::istream;
template <typename T>
istream& operator>>(istream& is, Polynomial<T> & p) {
    if (is.peek() == '[')
        is.get();
    else
        is.setstate(std::ios::failbit);
    typename Polynomial<T>::StorageT tempStorage;

    typename Polynomial<T>::CoefT cf;
    while ( is >> cf ) {
        tempStorage.push_back(cf);
    }

    p.setCoefs( tempStorage );
    return is;
}

template <typename T>
istream& operator>>(istream& is, Polynomial< Polynomial<T> > & p) {
*/

template<typename T>
void loadPolyFromString( Polynomial< Polynomial<T> > & p,
        std::string const & s ) {

    typedef Polynomial< Polynomial<T> > PolyT;

    static const boost::regex e("(\\[.*?\\] ?)*");

    boost::smatch mt;
    if (boost::regex_search(s, mt, e)) {

        typename PolyT::StorageT tempStorage;
        //std::copy(mt.captures.begin(), mt.end(), std::back_inserter(tempStorage));
        p.setCoefs( tempStorage );
    }

/*

    typename PolyT::CoefT cf;
    while ( is >> cf ) {
        tempStorage.push_back(cf);
    }

*/}

using std::ostream;
template <typename T>
ostream& operator<<(ostream& os, Polynomial<T> const & p) {
    os << '[';
    std::copy(p.getCoefs().begin(), p.getCoefs().end(),
            std::ostream_iterator<T>(os, " "));
    os << ']';
    return os;
}

/** Template type for multivariate polynomials
 (using recursive instantiation)
 @param VarCnt — polynomial variables count
 @param Coef — polynomial coefficient type */
template<int VarCnt, typename Coef>
class MVPolyType {
public:
    typedef Polynomial<typename MVPolyType<VarCnt - 1, Coef>::ResultT> ResultT;
};

/** Specialization of MVPolyType
 for stopping recursive instantiation */
template<typename Coef>
class MVPolyType<1, Coef> {
public:
    typedef Polynomial<Coef> ResultT;
};

using std::cout;
using std::endl;

int main() {
    typedef MVPolyType<2, int>::ResultT Poly2;
    Poly2 p2;
    std::string s("[[1 2 3][4 5 6]]");
    //std::istringstream iss();
    //std::istringstream iss("[4 5 6]");
    loadPolyFromString(p2, s);
    std::cout << "Число переменных Poly2: " << Poly2::VAR_CNT << std::endl;
    std::cout << "Тип Poly2: " << typeid(p2).name() << std::endl;
    std::cout << "p2: " << p2 << std::endl;

}
