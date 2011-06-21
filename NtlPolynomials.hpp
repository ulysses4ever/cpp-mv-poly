/*
 * NtlPolynomials.hpp
 *
 *  Created on: 21.06.2011
 *      Author: ulysses
 */

#ifndef NTLPOLYNOMIALS_HPP_
#define NTLPOLYNOMIALS_HPP_

#include <map>

#include <boost/foreach.hpp>

#include "mv_poly.hpp"
#include "Point.hpp"

/**
 * TODO: add comment
 */
template<template <typename> class OrderPolicy, typename T>
class PowerPrinter {
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
    PowerPrinter(CoefT const & x, Polynomial<T> const & p) :
            x(x) {
        data = polyToDegCoefMap<OrderPolicy>(p);
    }

    friend
    std::ostream & operator<<(
            std::ostream & os,
            PowerPrinter<OrderPolicy, T> const & pp) {
        //BOOST_FOREACH(typename PointCoefMap::value_type const & pt_cf, pp.data) {
        if (pp.data.empty())
            return os;
        for(typename PointCoefMap::const_iterator it = pp.data.begin();
                it != --pp.data.end(); ++it) {
            typename PointCoefMap::value_type const & pt_cf = *it;
            if (pt_cf.second != CoefT::zero()) {
                os << "a^" << pp.log(pt_cf.second)
                        << " X^" << pt_cf.first << " + ";
            }
        }
        typename PointCoefMap::value_type const & pt_cf = *(--pp.data.end());
        os << "a^" << pp.log(pt_cf.second)
                                << " X^" << pt_cf.first;
        return os;
    }
};

template<template <typename> class OrderPolicy, typename T>
PowerPrinter<OrderPolicy, T> makePowerPrinter(
        typename Polynomial<T>::CoefT const & x,
        Polynomial<T> const & p) {
    return PowerPrinter<OrderPolicy, T>(x, p);
}

#endif /* NTLPOLYNOMIALS_HPP_ */
