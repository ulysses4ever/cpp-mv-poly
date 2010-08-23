/** @file mv_poly.hpp
 *
 * Implementation of multivariate polynomials and accompanying routines.
 * The main idea of current implementation of multivariate polynomials
 * is in the fact that polynomial from say n variables is just the
 * polynomial from one variable with coefficients being (from the ring of)
 * polynomials from n-1 variables. This idea is implemented on the basis of
 * template recursive instantiation.
 *
 * @author Artem Pelenitsyn
 * @date 2010-06-28
 */
#ifndef MV_POLY_HPP_
#define MV_POLY_HPP_

#include <algorithm>
#include <functional>
#include <deque>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>

#include <tr1/functional>

#include "Point.hpp"
#include "CoefficientTraits.hpp"

namespace TestMVPoly {
void outputTest();
}

namespace mv_poly {

/// \cond
template<typename T>
class Polynomial;
/// \endcond

/** Neat template type for actually getting multivariate polynomials.
 * Delivers user from writing \c Polynomial<Polynomial<… Polynomial<int>… >.
 @param VarCnt — polynomial variables count
 @param Coef — polynomial coefficient type */
template<int VarCnt, typename Coef>
struct MVPolyType {
    typedef Polynomial<typename MVPolyType<VarCnt - 1, Coef>::ResultT> ResultT;
};

/// \cond
/* Specialization of MVPolyType
 for stopping recursive instantiation */
template<typename Coef>
struct MVPolyType<1, Coef> {
    // multivariate polynomial from 1 variable is just Polynomial
    typedef Polynomial<Coef> ResultT;
};
/// \endcond

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
    template<typename S>
    friend
    std::istream& operator>>(std::istream& is, Polynomial<S> & p);

    friend
    void TestMVPoly::outputTest();

    Polynomial() : data(1) /* no-argument constructor
                                    for data element means 0 */ {}

    explicit Polynomial(std::string const & s) {
        loadPolyFromString(*this, s );
    }

    /**
     * Template for computing the type of multivariate polynomial coefficients.
     */
    template<typename S>
    class PolyCoeff {
    public:
        typedef S Type;
    };

    template<typename S>
    class PolyCoeff< Polynomial<S> > {
    public:
        typedef typename PolyCoeff<S>::Type Type;
    };

    /**
     * Multivariate polynomial coefficients type.
     * For <tt>Polynomial<Polynomial<… Polynomial<T>… ></tt> it's \c T.
     */
    typedef typename PolyCoeff<T>::Type CoefT;

    /**
     * Number of polynomial variables. Actually it's nestedness depth of type
     * <tt>Polynomial<Polynomial<… Polynomial<T>… ></tt>.
     */
    static const int VAR_CNT = 1 + VarCnt<T>::result;

    /**
     * Subscript operator aimed for use with “plain” (one-variable) polynomials
     * (i.e.\ Polynomial<T> where T is not of the form Polynomial<S>).
     * @param idx index of (plain) polynomial coefficient to be returned.
     * @return (plain) polynomial coefficient with index \c idx.
     */
    CoefT operator[](int idx) const;

    /**
     * Subscript operator aimed for use with multivariate polynomials
     * (i.e.\ Polynomial<T> where T is not of the form Polynomial<S>, S is possibly
     * of the form Polynomial<U> and so on). Coefficients of mv-polynomial is
     * numbered with points (Point<Dim>), and point dimension (Dim) should be
     * equal to the number of polynomial variables.
     * @param pt point-index (sometimes called multi-indexd) of polynomial
     * coefficient to be returned.
     * @return polynomial coefficient with point-index \c pt.
     */
    template<template <typename PointImpl> class OrderPolicy>
    CoefT operator[](Point<VAR_CNT, OrderPolicy> const & pt) const;

    /**
     * Multiply polynomial on a scalar (assignment version).
     * @param c Scalar to multiply on.
     * @return This polynomial multiplyed on \c c.
     */
    Polynomial operator*=(CoefT const & c) {
        //using std::tr1::bind;
        //using std::tr1::placeholders::_1;
        //std::for_each(data.begin(), data.end(), // binding overloaded functions
                //bind(operator*=, _1, c)); // from different scopes is really messy
        for (typename StorageT::iterator it = data.begin(); it != data.end(); ++it) {
            (*it) *= c;
        }
        return *this;
    }

    /**
     * Multiply polynomial on a monomial represented by its degree (assignment version)
     * aimed for use with multivariate polynomials. If the polynomial has
     * VAR_CNT variables <tt>x = x_1 x_2 … x_{VAR_CNT}</tt>, then
     * <tt>Point<VAR_CNT> m</tt> represents monomial
     * <tt>x^m = x_1^{m_1} x_2^{m_2} … x_{VAR_CNT}^{m_{VAR_CNT}}</tt>
     * @param c Monomial x^{\c m} to multiply on.
     * @return This polynomial multiplyed on \c x^{\c m}.
     */
    template<template <typename PointImpl> class OrderPolicy>
    Polynomial operator<<=(Point<VAR_CNT, OrderPolicy> const & m);

    /**
     * Multiply polynomial on a monomial represented by its degree (assignment version)
     * aimed for use with “plain” (one-variable) polynomials.
     * @param m  x^{\c m} to multiply on.
     * @return This polynomial multiplyed on \c x^{\c m}.
     */
    Polynomial operator<<=(int m);

    /**
     * Polynomial addition (assignment version).
     * @param p[in] Polynomial to be added to this.
     * @return This polynomial after addition \c p.
     */
    Polynomial operator+=(Polynomial const & p);

    /**
     * Polynomial subtraction (assignment version).
     * @param p[in] Polynomial to be subtracted from this.
     * @return This polynomial after addition \c p.
     */
    Polynomial operator-=(Polynomial const & p);

    /**
     * Polynomial comparison for equality after normalization.
     * @param lhs Left-hand side operand for comparison
     * @param rhs Right-hand side operand for comparison
     * @return True iff all coefficients of both polynomials are equal, false
     * otherwise. If high coefficients of one polynomial are zero, they considered
     * equal to absent high coefficients of other polynomial.
     */
    template<typename S>
    friend
    inline
    bool operator==(Polynomial<S> const & lhs, Polynomial<S> const & rhs);

    /** As polynomial is actually a special container type, it has distinguished
     * type for it's elements.
     */
    typedef T                   ElemT;

    typedef std::deque<ElemT>   StorageT;

    StorageT const & getCoefs() const     { return data; }

    static Polynomial getId() {
        std::string strRep;
        std::fill_n(std::back_inserter(strRep), VAR_CNT, '[');
        strRep.push_back('1');
        std::fill_n(std::back_inserter(strRep), VAR_CNT, ']');
        return Polynomial(strRep);
    }

private:
    /**
     * Implementation detail: used in operator[](Point<VAR_CNT>) for going
     * deeper in the polynomial nest. If client calls
     * <tt>MVPolyType<n, T>::operator[Point<n> pt]</tt>,
     * then for all k we will use
     * <tt>MVPolyType<n - k, T>::operator[Slice<n, k>]</tt>,
     * where <tt>Slice<n, k></tt> offers access to the part of pt, namely
     * <tt>pt[k, k+1, ..., n-1]</tt>
     * @param sl
     * @return
     */
    template<int Dim, template <typename PointImpl> class OrderPolicy,
        int Offset>
    CoefT operator[](Slice<Dim, OrderPolicy, Offset> const & sl) const;

    template<typename T1, typename S, typename Pt>
    friend
    T1 applySubscript(S const & el, Pt const & pt);

    /**
     * Implementation detail: used in operator<<=(Point<VAR_CNT>) for going
     * deeper in the polynomial nest.
     * @param m Slice of the monomial got in initial client call for
     * <tt>operator<<=(Point<VAR_CNT>)</tt>.
     * @return This polynomial with one more dimension multiplied on the next
     * initial point component pointed to by the current slice.
     */
    template<int Dim, template <typename PointImpl> class OrderPolicy,
        int Offset>
    Polynomial operator<<=(Slice<Dim, OrderPolicy, Offset> const & m);

    template<typename S, typename Pt>
    friend
    void applyMonomialMultiplication(S & elem, Pt const & monomial);


    /**
     * Delete trailing zeros in coefficient collection \c data.
     */
    void normalization() {
        if (data.size() < 2)
            return; // we will not delete single zero
        ElemT tempDefElem = ElemT();
        typename StorageT::iterator it = --data.end(); // it is still valid
        do {
            if (*it == tempDefElem)
                --it;
            else
                break;
        } while (it != data.begin());
        if (it == data.begin()) // we'll not delete all the elements
            ++it;
        data.erase(it, data.end());
    }

    void setCoefs(StorageT const & data)  { this->data = data; }

    StorageT data;
};

template<typename T>
inline
bool operator==(Polynomial<T> const & lhs, Polynomial<T> const & rhs) {
       const_cast<Polynomial<T> &>(lhs).normalization();
       const_cast<Polynomial<T> &>(rhs).normalization();
       // yep, I know, that every time I use const_cast God kills a kitten
       return lhs.data == rhs.data;
}

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

/** We need generic subscript utility as with Polynomial::operator[] we
 * can't use template specialization (the rule for template members of
 * template class: specialize member only for enclosing class specialization
 * cf. [VJ, 12.3.3]).
 *
 * The point type here (Pt) is actually either Point or Slice, as opposed to
 * Polynomial::operator[] we don't handle the variants separately because it's
 * (some kind) hidden implementation detail — with Polynomial::operator[] we
 * wanna be as clear as we can, so parameteres types are as strict as they
 * could be.
 *
 * @param el the object to be subscripted
 * @param pt generic subscription index. The function could be tuned for speical
 * types of \c pt.
 * @return result of suscript operator applied to \c el at \c pt.
 */
template<typename T, typename S, typename Pt>
T applySubscript(S const & el, Pt const & pt) {
        return el[pt];
}

/// When slice is maximal Slice<Dim, Dim - 1> stop slicing a point.
template<typename T, typename S, int Dim,
    template <typename PointImpl> class OrderPolicy>
T applySubscript(S const & el, Slice<Dim, OrderPolicy, Dim - 1> const & pt) {
        return el[pt[0]];
}

/**
 * Handling Slice<Dim, Dim>. Actually this type is illegal
 * (in Slice<Dim, Offset> we should have Dim > Offset). If we get here, it means
 * client used MVPolyType<1, T>::operator[](Point<1>) instead of operator[](int).
 * But client is always right so we cope with this correctly.
 */
template<typename T, typename S, int Dim,
    template <typename PointImpl> class OrderPolicy>
T applySubscript(S const & el, Slice<Dim, OrderPolicy, Dim> const & pt) {
        return el;
}

template<typename T>
template<template <typename PointImpl> class OrderPolicy>
typename Polynomial<T>::CoefT
Polynomial<T>::operator[](Point<VAR_CNT, OrderPolicy> const & pt) const {
    if (pt[0] < (int)0 || (int)data.size() <= pt[0])
        return CoefT();
    else
        return applySubscript<Polynomial<T>::CoefT>(data[pt[0]], make_slice(pt));
}

template<typename T>
template<int Dim, template <typename PointImpl> class OrderPolicy, int Offset>
typename Polynomial<T>::CoefT
Polynomial<T>::operator[](Slice<Dim, OrderPolicy, Offset> const & sl) const {
    if (sl[0] < int(0) || int(data.size()) <= sl[0])
        return CoefT();
    else
        return applySubscript<Polynomial<T>::CoefT>(data[sl[0]], make_slice(sl));
}

template<typename T>
typename Polynomial<T>::CoefT
Polynomial<T>::operator[](int pt) const {
    if (pt < (int)0 || (int)data.size() <= pt)
        return CoefT();
    else
        return data[pt];
}

/**
 * Convolution-like operation for polynomials with immediate return of \c m\ -th
 * compoment (in general the result of convolution of two sequences is a
 * sequence also). We need a hint where to stop convoluting our two polynomial/
 * sequences — \c degf plays this role.
 * @param f One polynomial to convolute.
 * @param u Second polynomial to convolute.
 * @param degf The hint where to stop convolute.
 * @param m The index of resulting polynomial/sequence that we — only
 * virtually — get with our convolution operation.
 * @return
 */
template<typename T>
inline
typename Polynomial<T>::CoefT
conv(
        Polynomial<T> const & f,
        Polynomial<T> const & u,
        Point<Polynomial<T>::VAR_CNT> const & degf,
        Point<Polynomial<T>::VAR_CNT> const & m) {
    typename Polynomial<T>::CoefT res;
    Point<Polynomial<T>::VAR_CNT> i;
    while(i <= degf) {
        res += f[i] * u[i + m - degf];
        ++i;
    }
    return res;
}

/**
 * Input polynomial from stream, format: [[a b c][e f]].
 * @param[in,out] is Stream which contains data for creating polynomial.
 * @param[out] p Polynomial to store.
 * @return \c is (conventionally).
 */
template <typename T>
std::istream& operator>>(std::istream& is, Polynomial<T> & p) {
    char c = 0;
    is >> c;
    if ('[' != c) {
        is.setstate(std::ios::failbit);
        return is;
    }

    typename Polynomial<T>::StorageT tempStorage;
    typename Polynomial<T>::ElemT el;
    while ( is.peek() != ']' && is ) {
        is >> el;
        tempStorage.push_back(el);
    }

    if (is) {
        is.ignore(); // ignore trailing ']'
        if (!tempStorage.empty())
            p.setCoefs(tempStorage);
        else
            p = Polynomial<T>();
    }
    return is;
}

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

/**
 * Generic output for mv-polynomials, format: [[a b c][e f]].
 * @param[out] os Target output stream.
 * @param[in] p Polynomial to be outputed in \c os.
 * @return Output stream \c os after polynomial \c p have been
 * outputed to \c os (conventionally).
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, Polynomial<T> const & p) {
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

template<typename T>
std::string toString(Polynomial<T> const & p) {
    std::ostringstream os;
    os << p;
    return os.str();
}

/**
 * Multiply polynomial on a scalar.
 * @param c Scalar to multiply on.
 * @param p Polynomial to be multiplied.
 * @return This polynomial multiplyed on \c c.
 */
template<typename T>
inline
Polynomial<T> operator*(typename Polynomial<T>::CoefT const & c, Polynomial<T> p) {
    // p — possible pass-by-copy optimization
    return p *= c;
}

/**
 * Multiply polynomial on a scalar.
 * @param p Polynomial to be multiplied.
 * @param c Scalar to multiply on.
 * @return This polynomial multiplyed on \c c.
 */
template<typename T>
inline
Polynomial<T> operator*(Polynomial<T> p, typename Polynomial<T>::CoefT const & c) {
    return p *= c;
}

/**
 * The purpose of the function is essentially the same as with applySubscript:
 * to pull structure recursion on <tt>Polynomial<… Polynomial<T>… ></tt> type
 * deeper, in this case — applying operator<<= on the way.
 * @param elem
 * @param monomial
 */
template<typename T, typename Pt>
inline
void applyMonomialMultiplication(T & elem, Pt const & monomial) {
    elem <<= monomial;
}

template<typename T, int Dim, template <typename PointImpl> class OrderPolicy>
inline
void applyMonomialMultiplication(
        T & elem,
        Slice<Dim, OrderPolicy, Dim - 1> const & monomial) {
    elem <<= monomial[0];
}

template<typename T, int Dim, template <typename PointImpl> class OrderPolicy>
inline
void applyMonomialMultiplication(
        T & elem,
        Slice<Dim, OrderPolicy, Dim> const & monomial) {
    return;
}

template<typename T>
template<template <typename PointImpl> class OrderPolicy>
inline
Polynomial<T> Polynomial<T>::operator<<=(Point<VAR_CNT, OrderPolicy> const & monomial) {
    for (typename StorageT::iterator it = data.begin(); it != data.end(); ++it)
        applyMonomialMultiplication(*it, make_slice(monomial));
    std::fill_n(std::front_inserter(data), monomial[0], ElemT());
    return *this;
}

template<typename T>
template<int Dim, template <typename PointImpl> class OrderPolicy, int Offset>
inline
Polynomial<T>
Polynomial<T>::operator<<=(Slice<Dim, OrderPolicy, Offset> const & monomial) {
    using std::tr1::bind;
    using std::tr1::placeholders::_1;
    for (typename StorageT::iterator it = data.begin(); it != data.end(); ++it)
        applyMonomialMultiplication(*it, make_slice(monomial));
    std::fill_n(std::front_inserter(data), monomial[0], ElemT());
    return *this;
}

template<typename T>
inline
Polynomial<T> Polynomial<T>::operator<<=(int monomial) {
    std::fill_n(std::front_inserter(data), monomial, ElemT());
    return *this;
}


/**
 * Multiply polynomial on a monomial represented by its degree
 * aimed for use with multivariate polynomials. If the polynomial has
 * VAR_CNT variables <tt>x = x_1 x_2 … x_{VAR_CNT}</tt>, then
 * <tt>Point<VAR_CNT> m</tt> represents monomial
 * <tt>x^m = x_1^{m_1} x_2^{m_2} … x_{VAR_CNT}^{m_{VAR_CNT}}</tt>
 * @param p Polynomial to be multiplied.
 * @param c Monomial x^{\c m} to multiply on.
 * @return This polynomial multiplyed on \c x^{\c m}.
 */
template<typename T>
inline
Polynomial<T> operator<<(
        Polynomial<T> p,
        Point<Polynomial<T>::VAR_CNT> const & m) {
    return p <<= m;
}

/**
 * Multiply polynomial on a monomial represented by its degree
 * aimed for use with “plain” (one-variable) polynomials.
 * @param p Polynomial to be multiplied.
 * @param m  x^{\c m} to multiply on.
 * @return This polynomial multiplyed on \c x^{\c m}.
 */

template<typename T>
inline
Polynomial<T> operator<<(Polynomial<T> p, int m) {
    return p <<= m;
}

template<typename T>
Polynomial<T> Polynomial<T>::operator+=(Polynomial<T> const & p) {
    int degDiff = (this->data).size() - p.data.size();
    if (degDiff < 0) { // if p has greater (“one-dimensional”) degree then
        // first: copy tail of p to this
        typename StorageT::const_iterator commonPartDelimeter(p.data.begin());
        std::advance(commonPartDelimeter, data.size());
        std::copy(commonPartDelimeter, p.data.end(), std::back_inserter(data));
        // second: add the rest of p to this
        typename StorageT::iterator itThis = data.begin();
        for (
                typename StorageT::const_iterator itP = p.data.begin();
                itP != commonPartDelimeter;
                ++itThis, ++itP) {
            *itThis += *itP;
        }
    } else { // if this polynomial has greater degree then
        typename StorageT::iterator itThis = data.begin();
        for (
                typename StorageT::const_iterator itP = p.data.begin();
                itP != p.data.end();
                ++itThis, ++itP) {
            *itThis += *itP;
        }
    }
    return *this;
}

template<typename T>
inline
Polynomial<T> operator+(Polynomial<T> lhs, Polynomial<T> const & rhs) {
    return lhs += rhs;
}

template<typename T>
Polynomial<T> Polynomial<T>::operator-=(Polynomial<T> const & p) {
    return (*this +=
                CoefficientTraits<CoefT>::addInverse(
                        CoefficientTraits<CoefT>::multId())
                * p);
}

template<typename T>
inline
Polynomial<T> operator-(Polynomial<T> lhs, Polynomial<T> const & rhs) {
    return lhs -= rhs;
}

} // namespace mv_poly

#endif /* MV_POLY_HPP_ */
