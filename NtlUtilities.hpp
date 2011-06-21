/*
 * NtlUtilities.hpp
 *
 *  Created on: 21.06.2011
 *      Author: ulysses
 */

#ifndef NTLUTILITIES_HPP_
#define NTLUTILITIES_HPP_

#include <iostream>
#include <string>
#include <sstream>

#include <cmath>

#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

template<typename F>
struct NTLPrimeFieldTtraits {
    typedef NTL::ZZ_pE ExtField;
    typedef NTL::ZZ_pX PolyOver;

    // assuming p from ZZ_p fits int
    static int Char() {return NTL::to_int(F::modulus());}
};

template<>
struct NTLPrimeFieldTtraits<NTL::GF2> {
    typedef NTL::GF2E ExtField;
    typedef NTL::GF2X PolyOver;

    static int Char() {return 2;}
};


/**
 * Initialization of NTL extended field.
 * @param polyExtStr String representation of an irreducible polynomial over F
 * for extendent field construction construction.
 * @return Size of result (extended) field.
 */
template<typename F>
size_t initExtendedField(
        std::string const & polyExtStr) {
    typename NTLPrimeFieldTtraits<F>::PolyOver p;
    std::istringstream is(polyExtStr);
    is >> p;
    typedef typename NTLPrimeFieldTtraits<F>::ExtField EF;
    EF::init(p);
    return pow(NTLPrimeFieldTtraits<F>::Char(), deg(p));
}

/**
 * Builds primitive element of field ExtF providing that it was constructed
 * using _primitive_ polynomial.
 * @return primitive element of ExtF field.
 */
template<typename ExtF>
ExtF getPrimitive() {
    std::istringstream is("[0 1]");
    ExtF x;
    is >> x;
    return x;
}

/**
 *
 * @param x Field primitive element.
 * @param field_size Number of elements in the field.
 */
template<typename F>
void printFieldInPowers(F const & x, size_t field_size) {
    using std::cout;
    using std::endl;
    for (size_t i = 0; i < field_size - 1; ++i)
        cout << "a^" << i << " = " << power(x, i) << endl;
    cout << endl;
}


#endif /* NTLUTILITIES_HPP_ */
