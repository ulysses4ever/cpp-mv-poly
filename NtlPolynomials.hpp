/*
 * NtlPolynomials.hpp
 *
 *  Created on: 21.06.2011
 *      Author: ulysses
 */

#ifndef NTLPOLYNOMIALS_HPP_
#define NTLPOLYNOMIALS_HPP_

#include <map>
#include <string>

#include <NTL/GF2.h>
#include <NTL/ZZ_p.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ_pE.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>

#include "mv_poly.hpp"
#include "Point.hpp"

namespace mv_poly {

typedef boost::mpl::vector<NTL::GF2E, NTL::ZZ_pE> NtlExtFieldTypes;
typedef boost::mpl::vector<NTL::GF2, NTL::ZZ_p> NtlPrimeFieldTypes;

/**
 * Effector (cf. Eckel, TIC++ vol. 2, ch. 4) that takes a polynomial ‘p’
 * and a field primitive element ‘a’ to print ‘p’ in a “pretty” form:
 * a^k x^(m, n) + a^k' x^(m', n') + …
 */
template<template <typename> class OrderPolicy, typename T>
class PowerPolyPrinter {

    typedef typename Polynomial<T>::CoefT CoefT;

    typedef Point<Polynomial<T>::VAR_CNT, OrderPolicy> PointT;

    typedef std::map<PointT, CoefT> PointCoefMap;

    CoefT a;

    PointCoefMap data;

    int log(CoefT const & cf) const {
        int result = 0;
        CoefT pw = CoefficientTraits<CoefT>::multId();
        while (pw != cf) {
            mul(pw, pw, a);
            ++result;
        }
        return result;
    }

    template<typename CoefT>
    std::string
    coefToString(
            CoefT const & cf,
            typename boost::enable_if<
                    boost::mpl::contains<NtlExtFieldTypes, CoefT> >::type * = 0
            ) const {
        int l = log(cf);
        return 0 == l ? "" : "a^" + boost::lexical_cast<std::string>(l) + " ";
    }

    template<typename CoefT>
    std::string
    coefToString(
            CoefT const & cf,
            typename boost::enable_if<
                    boost::mpl::contains<NtlPrimeFieldTypes, CoefT> >::type * = 0
            ) const {
        return cf == CoefficientTraits<CoefT>::multId() ? "" :
                boost::lexical_cast<std::string>(cf) + " ";
    }

 public:

    PowerPolyPrinter(Polynomial<T> const & p, CoefT const & a) : a(a) {
        data = polyToDegCoefMap<OrderPolicy>(p);
    }

   friend
    std::ostream & operator<<(
            std::ostream & os,
            PowerPolyPrinter const & pp) {
        if (pp.data.empty())
            return os;
        for(typename PointCoefMap::const_iterator it = pp.data.begin();
                it != --pp.data.end(); ++it) {
            typename PointCoefMap::value_type const & pt_cf = *it;
            if (pt_cf.second != CoefT::zero()) {
                const std::string strCoef = pp.coefToString(pt_cf.second);
                if (pt_cf.first == PointT()) {
                    os << (strCoef == "" ? "1 " : strCoef) + "+ ";
                } else
                    os << strCoef
                            << "X^" << pt_cf.first << " + ";
            }
        }
        typename PointCoefMap::value_type const & pt_cf = *(--pp.data.end());
        os << pp.coefToString(pt_cf.second)
                                << "X^" << pt_cf.first;
        return os;
    }

};

template<template <typename> class OrderPolicy, typename T>
PowerPolyPrinter<OrderPolicy, T>
makePowerPrinter(
        Polynomial<T> const & p,
        typename Polynomial<T>::CoefT const & x
            = CoefficientTraits<typename Polynomial<T>::CoefT>::addId()) {
    return PowerPolyPrinter<OrderPolicy, T>(p, x);
}

} // namespace mv_poly

#endif /* NTLPOLYNOMIALS_HPP_ */
