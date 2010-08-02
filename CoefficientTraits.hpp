/**
 * @file CoefficientTraits.hpp
 *
 * Tuning set of polynomial coefficient-related operations.
 *
 * @author Artem Pelenitsyn
 * @date 2010-07-28
 */

#ifndef COEFFICIENTTRAITS_HPP_
#define COEFFICIENTTRAITS_HPP_

namespace mv_poly {

template<typename CoefT>
struct CoefficientTraits {
    static CoefT multInverse(CoefT const & c) {
        //return inv(c); // NTL-way // TODO: implement choosing inv for NTL types
        return c;
    }

    static CoefT addInverse(CoefT const & c) {
        //return inv(c); // NTL-way // TODO: implement choosing inv for NTL types
        return c;
    }

    static CoefT multId() {
        //return inv(c); // NTL-way // TODO: implement choosing inv for NTL types
        return CoefT();
    }
};

#include <NTL/GF2.h>
template<>
struct CoefficientTraits<NTL::GF2> {
    static NTL::GF2 multInverse(NTL::GF2 const & c) {
        return NTL::inv(c);
    }

    static NTL::GF2 addInverse(NTL::GF2 const & c) {
        return -c;
    }

    static NTL::GF2 multId() {
        return NTL::to_GF2(1);
    }
};

} // namespace mv_poly

#endif /* COEFFICIENTTRAITS_HPP_ */
