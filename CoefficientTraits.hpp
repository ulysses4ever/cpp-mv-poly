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

#include <boost/mpl/contains.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>

#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>

namespace mv_poly {

/**
 * \class CoefficientTraits
 * Provide basic algebraic properties of the
 * polynomial coefficients type.
 *
 * \note
 * We suppose that coefficient set forms a field,
 * so we want the type to have a way to obtain 0, 1, -a, a^{-1}.
 */
template<typename CoefT, typename Enable = void>
struct CoefficientTraits {

    /**
     * Obtaining multiplicative inverse in the field.
     * @param c Non-zero element.
     * @return \c c^{-1} (additive inverse in CoefT field).
     */
    static CoefT multInverse(CoefT const & c) {
        return 1 / c;
    }

    /**
     * Obtaining multiplicative inverse in the field.
     * @param c Any element of CoefT.
     * @return \c -c (additive inverse in CoefT field).
     */
    static CoefT addInverse(CoefT const & c) {
        return -c;
    }

    /**
     * Obtaining multiplicative identity of CoefT field.
     * @return Multiplicative identity (\c 1) of CoefT field.
     */
    static CoefT multId() {
        return 1;
    }

    /**
     * Obtaining additive identity of CoefT field.
     * @return Additive identity (\c 1) of CoefT field.
     */
    static CoefT addId() {
        return CoefT();
    }
};

typedef boost::mpl::vector<NTL::GF2E, NTL::ZZ_pE, NTL::GF2, NTL::ZZ_p>
    NtlFieldTypes;

/**
 * CoefficientTraits template specialization for NTL field types
 * (triggers with the boost::enable_if help).
 */
template<typename T>
struct CoefficientTraits<
        T,
        typename boost::enable_if<
            boost::mpl::contains<NtlFieldTypes, T>
        >::type > {

    static T multInverse(T const & c) {
        return NTL::inv(c);
    }

    static T addInverse(T const & c) {
        return -c;
    }

    static T multId() {
        T a;
        NTL::set(a);
        return a;
    }

    static T addId() {
        T zero;
        return zero;
    }
};

} // namespace mv_poly

#endif /* COEFFICIENTTRAITS_HPP_ */
