/** @file
 *
 * Implementation of multivariate polynomials and accompanying routines.
 * The main idea of current implementation of multivariate polynomials
 * is in the fact that polynomial from say n variables is just the
 * polynomial from one variable with coefficients being (from the ring of)
 * polynomials from n-1 variables. This idea is implemented on the basis of
 * template recursive instantiation.
 *
 * @author Artem Pelenitsyn
 */

#include <algorithm>
#include <functional>
#include <deque>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>

#include <tr1/array>

// #include <boost/logic/tribool.hpp>

/**
 * Point in N-dimensional integer lattice.
 * @param Dim Dimension of point lattice.
 */
template<int Dim>
class Point {
    typedef std::tr1::array<int, Dim> ImplType;

    ImplType data;

public:
    typename ImplType::reference
    operator[](typename ImplType::size_type n) {
        return data[n];
    }

    typename ImplType::const_reference
    operator[](typename ImplType::size_type n) const {
        return data[n];
    }

    typedef typename ImplType::iterator iterator;

    typedef typename ImplType::const_iterator const_iterator;

    iterator
    begin()
    { return data.begin(); }

    const_iterator
    begin() const
    { return data.begin(); }

    iterator
    end()
    { return data.end(); }

    const_iterator
    end() const
    { return data.end(); }

};

template<int Dim>
inline
bool byCoordinateLess(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return std::inner_product(lhs.begin(), rhs.begin(), true,
            std::logical_and<bool>(), std::less_equal<int>());
}

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
 * \note The class isn't meant to be used in application (client) code.
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
    // void loadPolyFromString(Polynomial<S> & p, std::string const & s );
    std::istream& operator>>(std::istream& is, Polynomial & p);

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

using std::istream;
/**
 * Input polynomial from stream;
 * @param[in,out] is Stream which contains data for creating polynomial.
 * @param[out] p Polynomial to store.
 * @return \c is (conventionally).
 */
template <typename T>
istream& operator>>(istream& is, Polynomial<T> & p) {
    char c = 0;
    is >> c;
    if ('[' != c) {
        is.setstate(std::ios::failbit);
        return is;
    }

    typename Polynomial<T>::StorageT tempStorage;
    typename Polynomial<T>::CoefT cf;
    while ( is.peek() != ']' && is ) {
        is >> cf;
        tempStorage.push_back(cf);
    }

    if (is) {
        is.ignore(); // ignore trailing ']'
        p.setCoefs( tempStorage );
    }
    return is;
}

/*
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
/**
 * Loading polynomial from the string.
 * @param[out] p Polynomial instance to get in loaded data.
 * @param[in] s String that defines contents of polynomial to de loaded.
 */
template<typename T>
void loadPolyFromString( Polynomial<T> & p, std::string const & s ) {
    std::istringstream is(s);
    is >> p;
}

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
