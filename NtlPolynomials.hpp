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
 * TODO: add comment
 */
template<template <typename> class OrderPolicy, typename T>
class PowerPolyPrinter {

    typedef typename Polynomial<T>::CoefT CoefT;

    typedef std::map<Point<Polynomial<T>::VAR_CNT, OrderPolicy>, CoefT>
        PointCoefMap;

    CoefT x;

    PointCoefMap data;

    int log(CoefT const & cf) const {
        int result = 0;
        CoefT pw = CoefficientTraits<CoefT>::multId();
        while (pw != cf) {
            mul(pw, pw, x);
            ++result;
        }
        return result;
    }

public:

    PowerPolyPrinter(Polynomial<T> const & p, CoefT const & x) : x(x) {
        data = polyToDegCoefMap<OrderPolicy>(p);
    }

    template<typename CoefT>
    std::string
    coef_tostr(
            CoefT const & cf,
            typename boost::enable_if<
                    boost::mpl::contains<NtlExtFieldTypes, CoefT> >::type * = 0
            ) const {
        return "a^" + boost::lexical_cast<std::string>(log(cf)) + " ";
    }

    template<typename CoefT>
    std::string
    coef_tostr(
            CoefT const & cf,
            typename boost::enable_if<
                    boost::mpl::contains<NtlPrimeFieldTypes, CoefT> >::type * = 0
            ) const {
        return cf == CoefficientTraits<CoefT>::multId() ? "" :
                boost::lexical_cast<std::string>(cf) + " ";
    }

    friend
    std::ostream & operator<<(
            std::ostream & os,
            PowerPolyPrinter<OrderPolicy, T> const & pp) {
        if (pp.data.empty())
            return os;
        for(typename PointCoefMap::const_iterator it = pp.data.begin();
                it != --pp.data.end(); ++it) {
            typename PointCoefMap::value_type const & pt_cf = *it;
            if (pt_cf.second != CoefT::zero()) {
                os << pp.coef_tostr(pt_cf.second)
                        << "X^" << pt_cf.first << " + ";
            }
        }
        typename PointCoefMap::value_type const & pt_cf = *(--pp.data.end());
        os << pp.coef_tostr(pt_cf.second)
                                << "X^" << pt_cf.first;
        return os;
    }

};

template<template <typename> class OrderPolicy, typename T>
PowerPolyPrinter<OrderPolicy, T> makePowerPrinter(
        Polynomial<T> const & p,
        typename Polynomial<T>::CoefT const & x
            = CoefficientTraits<typename Polynomial<T>::CoefT>::addId()) {
    return PowerPolyPrinter<OrderPolicy, T>(p, x);
}
} // namespace mv_poly

#endif /* NTLPOLYNOMIALS_HPP_ */
