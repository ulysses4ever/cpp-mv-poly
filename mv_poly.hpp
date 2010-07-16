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

#include "Point.hpp"


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
//    template<typename S>
//    friend
//    std::istream& operator>>(std::istream& is, Polynomial & p);
//
//    template<typename S>
//    friend
//    std::ostream& operator<<(std::ostream& os, Polynomial & p);

    Polynomial() {}

    Polynomial(std::string const & s) {
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

    typedef typename PolyCoeff<T>::Type CoefT;

    static const int VAR_CNT = 1 + VarCnt<T>::result;

    template<int Dim>
    CoefT operator[](Point<Dim> const & pt) const;

    template<int Dim, int Span>
    CoefT operator[](Slice<Dim, Span> const & sl) const;

    CoefT operator[](int pt) const;

    typedef T                   ElemT;

    typedef std::deque<ElemT>   StorageT;

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

template<typename T, typename S, typename Pt>
T apply_subscript(S const & el, Pt const & pt) {
        return el[pt];
}

template<typename T, typename S, int Dim>
T apply_subscript(S const & el, Slice<Dim, Dim - 1> const & pt) {
        return el[pt[0]];
}

//template<typename T>
//inline
//Polynomial<T>::CoefT
//apply_subscript(T const & el, Point<0> const & pt) {
//    return el;
//}

//CoefT operator[](Point<VAR_CNT> const & pt) const;
template<typename T>
template<int Dim>
typename Polynomial<T>::CoefT
Polynomial<T>::operator[](Point< Dim/*Polynomial<T>::VAR_CNT*/ > const & pt) const {
    //assert(Dim == VAR_CNT);
    if (pt[0] < (int)0 || (int)data.size() <= pt[0])
        return CoefT();
    else
        //return data[pt[0]][make_slice(pt)];
        return apply_subscript<Polynomial<T>::CoefT>(data[pt[0]], make_slice(pt));
}

template<typename T>
template<int Dim, int Span>
typename Polynomial<T>::CoefT
Polynomial<T>::operator[](Slice<Dim, Span> const & sl) const {
    if (sl[0] < 0 || data.size() <= sl[0])
        return CoefT();
    else
        return apply_subscript<Polynomial<T>::CoefT>(data[sl[0]], make_slice(sl));
}


template<typename T>
typename Polynomial<T>::CoefT
Polynomial<T>::operator[](int pt) const {
    if (pt < (int)0 || (int)data.size() <= pt)
        return CoefT();
    else
        return data[pt];
}

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
    typename Polynomial<T>::ElemT el;
    while ( is.peek() != ']' && is ) {
        is >> el;
        tempStorage.push_back(el);
    }

    if (is) {
        is.ignore(); // ignore trailing ']'
        p.setCoefs( tempStorage );
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

using std::ostream;

/**
 * Generic output for mv-polynomials.
 * @param[out] os Target output stream.
 * @param[in] p Polynomial to be outputed in \c os.
 * @return Output stream \c os after polynomial \c p have been
 * outputed to \c os (conventionally).
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

#endif /* MV_POLY_HPP_ */
