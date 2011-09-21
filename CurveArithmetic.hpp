/*
 * CurveArithmetic.hpp
 *
 *  Created on: 09.09.2011
 *      Author: ulysses
 */

#ifndef CURVEARITHMETIC_HPP_
#define CURVEARITHMETIC_HPP_

#include <vector>

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
    while (l-- > 0)
        result.push_back(p++);
    return result;
}

} // namespace mv_poly

#endif /* CURVEARITHMETIC_HPP_ */
