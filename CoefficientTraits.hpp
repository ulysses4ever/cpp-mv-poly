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

#include <loki/Typelist.h>
#include <loki/TypeManip.h>

#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>

namespace mv_poly {

template<typename CoefT>
struct DefaultCoefficientTraits {
    static CoefT multInverse(CoefT const & c) {
        return 1 / c;
    }

    static CoefT addInverse(CoefT const & c) {
        return -c;
    }

    static CoefT multId() {
        return 1;
    }
};

template<typename T>
struct NtlCoefficientTraits {
    static T multInverse(T const & c) {
        return NTL::inv(c);
    }

    static T addInverse(T const & c) {
        return -c;
    }

    static T multId() {
        static T a;
        clear(a);
        return a;
    }
};

template<typename T>
class isNtlType {
    typedef LOKI_TYPELIST_4(NTL::GF2, NTL::GF2E, NTL::ZZ_p, NTL::ZZ_pE)
            NtlFiniteFieldTypes;

public:
    enum { result = Loki::TL::IndexOf<NtlFiniteFieldTypes, T>::value >= 0 };
};

template<typename CoefT>
struct CoefficientTraits : public Loki::Select<
                            isNtlType<CoefT>::result,
                            NtlCoefficientTraits<CoefT>,
                            DefaultCoefficientTraits<CoefT> >::Result {};

} // namespace mv_poly

#endif /* COEFFICIENTTRAITS_HPP_ */
