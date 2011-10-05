/*
 * CurveArithmetic.hpp
 *
 *  Created on: 09.09.2011
 *      Author: ulysses
 */

#ifndef CURVEARITHMETIC_HPP_
#define CURVEARITHMETIC_HPP_

#include <array>
#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>

#include <iostream>

#include "Point.hpp"

namespace mv_poly {

// deadly want template aliases
// typedef std::vector< Point<N, OrderPolicy> > BasisCollection<N, OrderPolicy>;

template<int r>
std::vector< Point<2, WeightedOrder<r, r + 1>::template impl > >
getHermitianCodeBasis(int l) {
    typedef Point<2, WeightedOrder<r, r + 1>::template impl > Pt;
    std::vector< Pt > result;
    result.reserve(l);
    Pt p;
    std::generate_n(std::back_inserter(result), l,
            std::bind<Pt (Pt::*)(int)>((&Pt::operator++),
                    p,
                    42 /* dummy argument for operator++(int) */));
//    while (l-- > 0)
//        result.push_back(p++);
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

template<int r, typename CurvePoint, typename FieldElem>
std::vector<CurvePoint>
getPlainHermitianCurveRationalPoints() {
//    auto power = &FieldElemTraits<FieldElem>::power;
    FieldElem x = FieldElemTraits<FieldElem>::getPrimitive();
    std::vector<CurvePoint> result;
    FieldElem zero = FieldElemTraits<FieldElem>::addId();
    FieldElem id = FieldElemTraits<FieldElem>::multId();
    CurvePoint cp;

    cp[0] = zero;
    cp[1] = zero;
    result.push_back(cp);

    cp[0] = zero;
    cp[1] = id;
    do {
        if (isPlainHermitianCurveRationalPoint<FieldElem>(r, cp)) {
            result.push_back(cp);
        }
        cp[1] *= x;
    } while (cp[1] != id);

    cp[0] = id;
    do {
        cp[1] = id;
        do {
            if (isPlainHermitianCurveRationalPoint<FieldElem>(r, cp)) {
                result.push_back(cp);
            }
            cp[1] *= x;
        } while (cp[1] != id);

        cp[0] *= x;
    } while (cp[0] != id);

    return result;
}

/**
 * Computes p^m ( = p_1^m_1 * p_2^m_2 * ... * p_n^m_n).
 */
template<typename FieldElem, typename Monom, typename CurvePoint>
FieldElem computeMonomAtPoint(Monom const & m, CurvePoint const & p) {

    return std::inner_product(p.begin(), p.end(), m.begin(),
            FieldElemTraits<FieldElem>::multId(),
            std::multiplies<FieldElem>(),
            FieldElemTraits<FieldElem>::template power<typename Monom::value_type>);
}

template<int r, typename FieldElem>
struct HermitianCodeTraits {

    //template<typename FieldElem, typename CurvePoint>
    struct impl {
        static auto getCodeBasis(int l) -> decltype(getHermitianCodeBasis<r>(l)) {
            return getHermitianCodeBasis<r>(l);
        }

        typedef decltype(getCodeBasis(42)) BasisCollection;

        typedef typename BasisCollection::value_type BasisElem;

        typedef std::array<FieldElem, 2> CurvePoint;

//        template<typename CurvePoint, typename FieldElem>
        static
        std::vector<CurvePoint>
        getRationalPoints() {
            return getPlainHermitianCurveRationalPoints<r, CurvePoint, FieldElem>();
        }

    };
};

} // namespace mv_poly

#endif /* CURVEARITHMETIC_HPP_ */
