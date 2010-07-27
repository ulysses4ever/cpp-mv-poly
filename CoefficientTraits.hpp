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

template<typename CoefT>
struct CoefficientTraits {
    CoefT inverse(CoefT const & c) {
        return inv(c); // NTL-way
    }
};

#endif /* COEFFICIENTTRAITS_HPP_ */
