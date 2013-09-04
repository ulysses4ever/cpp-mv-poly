/*
 * NtlUtilities.hpp
 *
 *  Created on: 21.06.2011
 *      Author: ulysses
 */

#ifndef NTLUTILITIES_HPP_
#define NTLUTILITIES_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <stdexcept>

#include <cmath>

#include <boost/lexical_cast.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>

#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>


typedef boost::mpl::vector<NTL::GF2E, NTL::ZZ_pE, NTL::GF2, NTL::ZZ_p>
    NtlFieldTypes;


// to be refactored out into separate file
template<typename FieldElem, typename Enable = void>
struct FieldElemTraits {
    template <typename DecimalType>
    static FieldElem power(FieldElem const & a, DecimalType n) {
        throw std::logic_error("NIY");
    }

    static FieldElem multId() {
        throw std::logic_error("NIY");
    }

    static FieldElem addId() {
        throw std::logic_error("NIY");
    }
};

template<typename FieldElem>
struct FieldElemTraits<
        FieldElem,
        typename boost::enable_if<
            boost::mpl::contains<NtlFieldTypes, FieldElem>
        >::type > {

    template <typename DecimalType>
    static FieldElem power(FieldElem const & a, DecimalType n) {
        return NTL::power(a, n);
    }

    static FieldElem multId() {
        FieldElem a;
        NTL::set(a);
        return a;
    }

    static FieldElem addId() {
        FieldElem zero;
        return zero;
    }

    static FieldElem getPrimitive() {
        if (!prim_ptr.get())
            throw std::logic_error("Primitive element is not initialized!");
        return *prim_ptr;
    }

    static void setPrimitive(FieldElem const & p) {
        prim_ptr.reset(new FieldElem(p));
//        isPrimitiveInitialized = true;
    }

private:
//    static bool isPrimitiveInitialized;

    static std::shared_ptr<FieldElem> prim_ptr;

};

//template<typename FieldElem>
//bool
//FieldElemTraits<
//        FieldElem,
//        typename boost::enable_if<
//            boost::mpl::contains<NtlFieldTypes, FieldElem>
//        >::type >::isPrimitiveInitialized = false;

template<typename FieldElem>
std::shared_ptr<FieldElem>
FieldElemTraits<
        FieldElem,
        typename boost::enable_if<
            boost::mpl::contains<NtlFieldTypes, FieldElem>
        >::type >::prim_ptr;

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
    FieldElemTraits<EF>::setPrimitive(getPrimitive<EF>());
    return pow(NTLPrimeFieldTtraits<F>::Char(), deg(p));
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

typedef boost::mpl::vector<NTL::GF2E, NTL::ZZ_pE> NtlExtFieldTypes;
typedef boost::mpl::vector<NTL::GF2, NTL::ZZ_p> NtlPrimeFieldTypes;

/**
 * Print NTL field elements as powers of a given primitive element.
 */
template<typename T>
class NtlPowerPrinter {

    typedef T ElemT;

    const ElemT a;

    const ElemT elem;

    int log(ElemT const & cf) const {
        int result = 0;
//        using std::cout;
//        using std::endl;
//        cout << "log_{ " << a << " }( " << cf << " )";
        ElemT pw = FieldElemTraits<ElemT>::multId();
        while (pw != cf) {
//            cout << " pw: " << pw;
//            cout << " cf: " << cf;
            mul(pw, pw, a);
            ++result;
        }
//        cout << "pw: " << pw;
        return result;
    }

    template<typename CoefT>
    std::string
    elemToString(
            CoefT const & cf,
            typename boost::enable_if<
                    boost::mpl::contains<NtlExtFieldTypes, CoefT> >::type * = 0
            ) const {
        if (cf == FieldElemTraits<CoefT>::addId())
            return "0";
        int l = log(cf);
        return 0 == l ? "1" : "a^" + boost::lexical_cast<std::string>(l);
    }

    template<typename CoefT>
    std::string
    elemToString(
            CoefT const & cf,
            typename boost::enable_if<
                    boost::mpl::contains<NtlPrimeFieldTypes, CoefT> >::type * = 0
            ) const {
        return boost::lexical_cast<std::string>(cf);
    }

 public:

    NtlPowerPrinter(ElemT const & elem) :
        a(FieldElemTraits<ElemT>::getPrimitive()), elem(elem) {}

    friend
    std::ostream & operator<<(
            std::ostream & os,
            NtlPowerPrinter const & pp) {
        os << pp.elemToString(pp.elem);
        return os;
    }

}; // class NtlPowerPrinter

template<typename ElemT>
inline
NtlPowerPrinter<ElemT> makeNtlPowerPrinter(ElemT const & elem) {
    return NtlPowerPrinter<ElemT>(elem);
}

#endif /* NTLUTILITIES_HPP_ */
