/*
 * CurveArithmetic.hpp
 *
 *  Created on: 09.09.2011
 *      Author: ulysses
 */

#ifndef CURVEARITHMETIC_HPP_
#define CURVEARITHMETIC_HPP_

#include <functional>
#include <iterator>
#include <vector>

#include <iostream>

#include "Point.hpp"

namespace mv_poly {

// deadly want template aliases
// typedef std::vector< Point<N, OrderPolicy> > BasisCollection<N, OrderPolicy>;

template<int r> //, typename OrderPolicy
std::vector< Point<2, WeightedOrder<r, r + 1>::template impl > >
getBasis(int l) {
    typedef Point<2, WeightedOrder<r, r + 1>::template impl > Pt;
    std::vector< Pt > result;
    result.reserve(l);
    Pt p;
   // std::function<Pt(int)> f(&Pt::operator++);
    //std::generate_n(std::back_inserter(result), l,
//            std::bind<Pt, Pt, int>(&Pt::operator++, p);//);
//            std::bind(std::plus<int>(), 3);
    while (l-- > 0)
        result.push_back(p++);
    return result;
}

template<typename FieldElem, typename CurvePoint>
bool isPlainHermitianCurveRationalPoint(int r, CurvePoint cp) {
    // Hermitian curve equation --v
    return FieldElemTraits<FieldElem>::power(cp[0], r + 1) -
            FieldElemTraits<FieldElem>::power(cp[1], r) - cp[1] ==
                    FieldElemTraits<FieldElem>::addId();
    // TODO: repetitive computations of x^n are ineffective
}

template<typename CurvePoint, typename FieldElem>
std::vector<CurvePoint>
getPlainHermitianCurveRationalPoints(int r, FieldElem x) {
//    auto power = &FieldElemTraits<FieldElem>::power;
    std::vector<CurvePoint> result;
    FieldElem zero = FieldElemTraits<FieldElem>::addId();
    CurvePoint cp;
    cp[0] = x;
    do {
        cp[1] = x;
        do {
            if (isPlainHermitianCurveRationalPoint<FieldElem>(r, cp)) {
                result.push_back(cp);
            }
            cp[1] *= x;
        } while (cp[1] != x);

        cp[0] *= x;
    } while (cp[0] != x);

    cp[0] = zero;
    cp[1] = x;
    do {
        if (isPlainHermitianCurveRationalPoint<FieldElem>(r, cp)) {
            result.push_back(cp);
        }
        cp[1] *= x;
    } while (cp[1] != x);

    cp[0] = zero;
    cp[1] = zero;
    result.push_back(cp);
    return result;
}

} // namespace mv_poly

#endif /* CURVEARITHMETIC_HPP_ */
