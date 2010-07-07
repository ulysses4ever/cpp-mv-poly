/** @file
 * Implementation of multivariate polynomials and accompanying routines.
 * @author Artem Pelenitsyn
 * The main idea of current implementation of multivariate polynomials
 * is in the fact that polynomial from say n variables is just the
 * polynomial from one variable with coefficients being (from the ring of)
 * polynomials from n-1 variables.
 */

#include <algorithm>
#include <deque>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <typeinfo>

#include <boost/regex.hpp>

/**
 * Variables counter for Polynomial class template.
 * Uses simple temp
 * late recursion. The actual formula for the number of
 * variables in Polynomial<T> is:
 * <tt>VarCnt< Polynomial<T> > = 1 + VarCnt<T> ></tt>
 * with VarCnt<T> is equal to 0 for all T not in form Polynomial<S>.
 *
 * For the purpose of counting the recurison depth in recursive template types
 * like Polynomial Vandevoorde&Josuttis suggest making enum constant
 * in recursive type itself, but we notice
 * some code bloating when defining template specialization for
 * recursion stopping. So we decided to make separate class (VarCnt) for this.
 *
 * \note The class isn't meant to be used in application code.
 */
template<typename T>
class VarCnt {
public:
    static const int result = 0;
};

/**
 * Generic polynomial class template.
 */
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

/// \cond
/*
 * Template recursion for defining polynomial variables counter continued.
 */
template<typename T>
class VarCnt< Polynomial<T> > {
public:
    static const int result = Polynomial<T>::VAR_CNT;
};
/// \endcond

/** Neat template type for actually getting multivariate polynomials.
 * Delivers user from writing \c Polynomial<Polynomial< ... Polynomial<int>...>.
 @param VarCnt — polynomial variables count
 @param Coef — polynomial coefficient type */
template<int VarCnt, typename Coef>
class MVPolyType {
public:
    typedef Polynomial<typename MVPolyType<VarCnt - 1, Coef>::ResultT> ResultT;
};

/// \cond
/* Specialization of MVPolyType
 for stopping recursive instantiation */
template<typename Coef>
class MVPolyType<1, Coef> {
public:
    // multivariate polynomial from 1 variable is just Polynomial
    typedef Polynomial<Coef> ResultT;
};
/// \endcond

/**
 * Loading polynomial from the string.
 * @param[out] p Polynomial instance to get in loaded data.
 * @param[in] s String that defines contents of polynomial to de loaded.
 */
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
*/

/*

    typename PolyT::CoefT cf;
    while ( is >> cf ) {
        tempStorage.push_back(cf);
    }
}

*/

using std::ostream;

/**
 * Generic output for mv-polynomials.
 * @param[out] os Target output stream.
 * @param[in] p Polynomial to be outputed in \c os.
 * @return Output stream \c os after polynomial \c p have been
 * outputed to \c os.
 */
template <typename T>
ostream& operator<<(ostream& os, Polynomial<T> const & p) {
    if (p.getCoefs().empty()) {
        os << "[]";
        return os;
    }
    os << '[';
    std::copy(p.getCoefs().begin(), --(p.getCoefs().end()),
            std::ostream_iterator<T>(os, " "));
    os << *(--(p.getCoefs().end()));
    os << ']';
    return os;
}
