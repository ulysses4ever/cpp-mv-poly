#include <algorithm>
#include <list>
#include <map>
#include <iterator>
#include <string>
#include <sstream>

#include <tr1/array>

#include <NTL/ZZ_p.h>
#include <NTL/GF2.h>

#define USE_TR1

#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"

#include "mv_poly.hpp"
#include "Point.hpp"
#include "bmsa.hpp"
#include "bmsa-decoding.hpp"
#include "NtlUtilities.hpp"
#include "NtlPolynomials.hpp"
#include "CurveArithmetic.hpp"

namespace TestMVPoly {

using namespace mv_poly;

using std::string;
using std::ostringstream;

void ioTestFor1Poly() {
    typedef MVPolyType<1, int>::ResultT Poly1;
    Poly1 p1;
    string s("[1 2 3]");
    ostringstream os;
    loadPolyFromString(p1, s);
    os << p1;
    ASSERT_EQUAL(s, os.str());

    Poly1 emptyPoly;
    s = "[]";
    os.str("");
    loadPolyFromString(emptyPoly, s);
    os << emptyPoly;
    ASSERT_EQUAL("[0]", os.str());
}

void outputTest() {
    typedef MVPolyType<1, int>::ResultT Poly1;
    typedef MVPolyType<2, int>::ResultT Poly2;

    // prepare Poly2 instance to output
    Poly1 p11("[1 2 3]"), p12("[3 2 1]"), p13("[1]");
    Poly2::StorageT st2;
    st2.push_back(p11);
    st2.push_back(p12);
    st2.push_back(p13);
    Poly2 p2;
    p2.setCoefs(st2);

    ostringstream os;
    os << p2;
    string s("[[1 2 3] [3 2 1] [1]]");
    ASSERT_EQUAL(s, os.str());
}

void inputTestForNPoly() {
    MVPolyType<2, int>::ResultT p2;
    loadPolyFromString(p2, "[[1 2 3] [3 2 1] [1]]");
    ostringstream os;
    os << p2;
    string s("[[1 2 3] [3 2 1] [1]]");
    ASSERT_EQUAL(s, os.str());

    MVPolyType<3, int>::ResultT p3;
    loadPolyFromString(p3, "[[[1 2] [3]] [[3] [2 1]] [[1]]]");
    os.str("");
    os << p3;
    s = "[[[1 2] [3]] [[3] [2 1]] [[1]]]";
    ASSERT_EQUAL(s, os.str());
}

void inputTestForNPolyOverGF() {
    NTL::ZZ_p::init(NTL::to_ZZ(2));
    MVPolyType<3, NTL::ZZ_p>::ResultT p3;
    loadPolyFromString(p3, "[[[1 2] [3]] [[3] [2 1]] [[1]]]");
    ostringstream os;
    os.str("");
    os << p3;
    string s = "[[[1 0] [1]] [[1] [0 1]] [[1]]]";
    ASSERT_EQUAL(s, os.str());
}

void polySubscript() {
    MVPolyType<2, int>::ResultT p;
    loadPolyFromString(p, "[[0 1 0] [1 0] [0] [1]]");
    Point<2> i, pt;
    pt[0] = 2; pt[1] = 1;
    ostringstream os;
    for ( ; i < pt; ++i) {
        os << p[i] << " ";
    }
    ASSERT_EQUAL("0 1 1 0 0 0 1 ", os.str());
}

void pointComparison() {
    Point<3> pt1, pt2, pt3;
    pt1[0] = 3; pt1[1] = 1; pt1[2] = 2;
    pt2[0] = 2; pt2[1] = 1; pt2[2] = 0;
    pt3[0] = 2; pt3[1] = 2; pt3[2] = 0;
    ASSERT(byCoordinateLess(pt2, pt1));  // pt1 <_p pt2
    ASSERT(!byCoordinateLess(pt1, pt2)); // pt2 \not <_p pt1 as pt1 <_p pt2 (see above)
    ASSERT(!byCoordinateLess(pt1, pt3)); // pt1 and pt3 incomparable

    Point<2> pt4, pt5;
    pt4[0] = 1; pt4[1] = 0;
    pt5[0] = 0; pt5[1] = 1;
    ASSERT(pt4 < pt5); // less by antilex
    ASSERT(!(pt5 < pt4));
    pt5[0] = 2; pt5[1] = 0;
    ASSERT(pt4 < pt5); // less by grading
    ASSERT(!(pt5 < pt4));

    ASSERT(!(pt5 < pt5));

    pt4[0] = 0; pt4[1] = 1;
    ASSERT(pt4 < pt5);
    ASSERT(!(pt5 < pt4));
}

void pointIncreasing() {
    Point<3> pt;
    pt[0] = 0; pt[1] = 0; pt[2] = 0;
    ++pt;
    ASSERT(pt[0] == 1 && pt[1] == 0 && pt[2] == 0);
    ++pt;
    ASSERT(pt[0] == 0 && pt[1] == 1 && pt[2] == 0);
    ++pt;
    ASSERT(pt[0] == 0 && pt[1] == 0 && pt[2] == 1);
    ++pt;
    ASSERT(pt[0] == 2 && pt[1] == 0 && pt[2] == 0);
    ++pt;
    ASSERT(pt[0] == 1 && pt[1] == 1 && pt[2] == 0);
    ++pt;
    ASSERT(pt[0] == 0 && pt[1] == 2 && pt[2] == 0);
    ++pt;
    ASSERT(pt[0] == 1 && pt[1] == 0 && pt[2] == 1);
}

void pointCollectionOperations() {
    Point<2> pt;
    std::list<Point<2> > s, sn, sig;
    pt[0] = 0; pt[1] = 1;
    s.push_back(pt);
    pt[0] = 2; pt[1] = 0;
    s.push_back(pt);
    pt[0] = 1; pt[1] = 0;
    s.push_back(pt);
    sn = getPartialMaximums(s);

    ASSERT_EQUAL(2, sn.size());

    pt[0] = 0; pt[1] = 1;
    ASSERT(std::find(sn.begin(), sn.end(), pt) != sn.end());

    pt[0] = 2; pt[1] = 0;
    ASSERT(std::find(sn.begin(), sn.end(), pt) != sn.end());

    sig = getConjugatePointCollection(sn);

    ASSERT_EQUAL(3, sig.size());

    pt[0] = 3; pt[1] = 0;
    ASSERT(std::find(sig.begin(), sig.end(), pt) != sig.end());

    pt[0] = 1; pt[1] = 1;
    ASSERT(std::find(sig.begin(), sig.end(), pt) != sig.end());

    pt[0] = 0; pt[1] = 2;
    ASSERT(std::find(sig.begin(), sig.end(), pt) != sig.end());

    pt[0] = 1; pt[1] = 1;
    s.push_back(pt);
    sn = getPartialMaximums(s);

    ASSERT_EQUAL(sn.size(), 2);

    pt[0] = 1; pt[1] = 1;
    ASSERT(std::find(sn.begin(), sn.end(), pt) != sn.end());

    pt[0] = 2; pt[1] = 0;
    ASSERT(std::find(sn.begin(), sn.end(), pt) != sn.end());

    sig = getConjugatePointCollection(sn);

    ASSERT_EQUAL(3, sig.size());

    pt[0] = 3; pt[1] = 0;
    ASSERT(std::find(sig.begin(), sig.end(), pt) != sig.end());

    pt[0] = 2; pt[1] = 1;
    ASSERT(std::find(sig.begin(), sig.end(), pt) != sig.end());

    pt[0] = 0; pt[1] = 2;
    ASSERT(std::find(sig.begin(), sig.end(), pt) != sig.end());
}

// follows Sakata [S88] example
void convolutionTest() {
    MVPolyType<2, NTL::GF2>::ResultT f, u;
    loadPolyFromString(u, "[[0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]");
    loadPolyFromString(f, "[[1 1] [1]]");
    Point<2> degf;
    degf[0] = 0;
    degf[1] = 1;
    Point<2> m;
    m[0] = 0;
    m[1] = 2;
    ASSERT_EQUAL(0, conv(f, u, degf, m));
    m[0] = 2;
    m[1] = 1;
    ASSERT_EQUAL(1, conv(f, u, degf, m));
}

void scalarMultiplication() {
    MVPolyType<2, int>::ResultT p("[[1 0 1] [1 1]]");
    p *= 2;
    ostringstream os;
    os << p;
    ASSERT_EQUAL("[[2 0 2] [2 2]]", os.str());
}

void monomialMultiplication() {
    ostringstream os;
    string init("[[1 0 1] [1 1]]");
    MVPolyType<2, int>::ResultT p(init);
    Point<2> pt;
    pt[0] = 0; pt[1] = 1;
    os << (p << pt);
    ASSERT_EQUAL("[[0 1 0 1] [0 1 1]]", os.str());
    os.str("");

    pt[0] = 0; pt[1] = 0;
    ASSERT_EQUAL(init, toString(p << pt));

    pt[0] = 1; pt[1] = 0;
    os << (p << pt);
    ASSERT_EQUAL("[[0] [1 0 1] [1 1]]", os.str());
}

void summation() {
    MVPolyType<2, int>::ResultT p("[[1 0 1] [1 1]]"), q("[[2 3] [0 2] [3]]");
    ASSERT_EQUAL(2*p, p + p);
    ASSERT_EQUAL("[[3 3 1] [1 3] [3]]", toString(p + q));
}

void equality() {
    MVPolyType<2, int>::ResultT p, q("[[0 0] [0]]");
    ASSERT_EQUAL(p, q);
}

void eval() {
    mv_poly::MVPolyType<1, int>::type p("[3 2 1]");
    ASSERT_EQUAL(11, p(2));

    mv_poly::MVPolyType<2, int>::type p1("[[3 2] [1]]");
    std::vector<int> pt(2);
    pt[0] = 1; pt[1] = 2;
    ASSERT_EQUAL(8, p1(pt));
}

void sakatasExample2D() {
    ostringstream os;
    typedef MVPolyType<2, NTL::GF2>::ResultT PolyT;
    PolyT u("[[0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]");
    Point<2> pt;
    pt[0] = 4; pt[1] = 1;
    BMSAlgorithm< PolyT > alg(u, pt);
    BMSAlgorithm< PolyT >::PolynomialCollection minset2 = alg.computeMinimalSet();
    copy(minset2.begin(), minset2.end(), std::ostream_iterator<PolyT>(os, "\n"));
    ASSERT_EQUAL(
            "[[1 0] [1 1] [0]]\n"
            "[[1 0 1] [1 1] [1]]\n"
            "[[1 1] [1 0] [0] [1]]\n", os.str());
}

void sakatasExample3D() {
    ostringstream os;
    Point<3> ptt;
    ptt[0] = 5; ptt[1] = 0; ptt[2] = 1;
    typedef MVPolyType<3, NTL::GF2>::ResultT PolyT3;
    PolyT3 v(
            "[[[1 1 1 1 0 0] [0 1 0 1 0] [1 1 0 0] [0 1 0] [0 0] [0] [1]]"//x = 0
            "[[1 1 0 1 1] [1 0 1 1] [0 1 1] [1 1] [1] [0]]" // x = 1
            "[[0 1 0 0] [0 0 1] [0 0] [1] [0]]" // x = 2
            "[[1 1 0] [1 0] [0] [1]] [[1 1] [0] [1]] [[1] [1]] [[0]]]" // x = 3-6
            );
    BMSAlgorithm< PolyT3 > alg3(v, ptt);
    BMSAlgorithm< PolyT3 >::PolynomialCollection minset = alg3.computeMinimalSet();

    os.str("");
    copy(minset.begin(), minset.end(), std::ostream_iterator<
            BMSAlgorithm< PolyT3 >::PolynomialCollection::value_type>(os, "\n"));

    ASSERT_EQUAL(
            "[[[1 1] [1]] [[0]] [[1]]]\n"
            "[[[0 1] [0 1] [0]] [[0 0] [0]] [[1]]]\n"
            "[[[1 1 1] [1] [1]] [[0 0] [0]] [[1]]]\n"
            "[[[1 0] [0 0] [1] [1]] [[0 0] [0] [0]] [[1] [1]] [[0]]]\n",
            os.str());
}

void testPolyToDegCoefMapConversion() {
    using namespace std;
    Polynomial<int> p;
    istringstream iss("[3 2 3]");
    iss >> p;
    auto res = polyToDegCoefMap< GradedAntilexMonomialOrder >(p);
    typedef decltype(res) DegCoefMap;
    ostringstream os;
    copy(res.begin(), res.end(), ostream_iterator<DegCoefMap::value_type>(os, ", "));
    ASSERT_EQUAL("(0, 3), (1, 2), (2, 3), ", os.str());

    MVPolyType<2, int>::type p1;
    iss.str("[[3 2] [3 1] [1]]");
    iss >> p1;
    auto res1 = polyToDegCoefMap< GradedAntilexMonomialOrder >(p1);
    typedef decltype(res1) DegCoefMap1;
    os.str("");
    copy(res1.begin(), res1.end(),
            ostream_iterator<DegCoefMap1::value_type>(os, ", "));
    ASSERT_EQUAL(
            "((0, 0), 3), ((1, 0), 3), ((0, 1), 2), ((2, 0), 1), ((1, 1), 1), ",
            os.str());
}

void polyPowerPrinting() {
    typedef NTL::GF2 PrimeField;
    typedef typename NTLPrimeFieldTtraits<PrimeField>::ExtField ExtField;

    std::ostringstream os;

    initExtendedField<PrimeField>("[1 1 1]");
    ExtField x = getPrimitive<ExtField>();

    MVPolyType<2, ExtField>::type
        p("[[[1]] [[1 1]]]"); //[1 0 1 1] [0 1 1] [1 1] [1] [0]

    os << makePowerPrinter< GradedAntilexMonomialOrder >(p);
    ASSERT_EQUAL("1 + a^2 X^(1, 0)", os.str());

    os.str("");
    MVPolyType<2, PrimeField>::type p1("[[1 1] [1 1 1]]");

    os << makePowerPrinter< GradedAntilexMonomialOrder >(p1);
    ASSERT_EQUAL("1 + X^(1, 0) + X^(0, 1) + X^(1, 1) + X^(1, 2)",
            os.str());
}

void curveArithmetic() {
    typedef NTL::GF2 PrimeField;
    typedef typename NTLPrimeFieldTtraits<PrimeField>::ExtField ExtField;

    initExtendedField<PrimeField>("[1 1 0 0 1]");
    ExtField x = FieldElemTraits<ExtField>::getPrimitive();


    auto cpts = getPlainHermitianCurveRationalPoints
            < 4, std::tr1::array<ExtField, 2>, ExtField
            >();
    typedef decltype(cpts) PointsCont;
    typedef PointsCont::value_type CPt;
    ASSERT_EQUAL(cpts.size(), 64);

//    cout << makeNtlPowerPrinter(x, cpts[0][0]) << ", "//<< cpts[0][1];
//            << makeNtlPowerPrinter(x, cpts[0][1]) << endl;
//
    auto b = getHermitianCodeBasis<4>(16);
    ASSERT_EQUAL(b.size(), 16);
    ASSERT(!(b.back()[0] == 0 && b.back()[1] == 0));
//
////    cout << b[0][0] << ", " << b[0][1] << endl;
//    cout << b[1][0] << ", " << b[1][1] << endl;
////    cout << b[2][0] << ", " << b[2][1] << endl;
//
//    cout << "computeMonomAtPoint(): " << makeNtlPowerPrinter(x,
//            computeMonomAtPoint<ExtField>(b[1], cpts[0])) << endl;
//
//    cout << makeNtlPowerPrinter(x,
//            FieldElemTraits<ExtField>::power<long>(cpts[0][0], b[1][0])) << endl;
//    cout << //makeNtlPowerPrinter(x,
//            FieldElemTraits<ExtField>::power<long>(cpts[0][1], b[1][1])/*)*/ << endl;
    std::tr1::array<int, 2> m;
    m[0] = 4; m[1] = 0;
    CPt cp;
    cp[0] = x;
    cp[1] = ExtField();
    ASSERT_EQUAL(computeMonomAtPoint<ExtField>(m, cp), // x^4 = x + 1
            NTL::power(x, 1) + NTL::power(x, 0));
}

void bmsaDecodingCLOS05Example() {
    typedef NTL::GF2 PrimeField;
    typedef typename NTLPrimeFieldTtraits<PrimeField>::ExtField ExtField;

    initExtendedField<PrimeField>("[1 1 1]");
    //ExtField a = FieldElemTraits<ExtField>::getPrimitive();
    FieldElemTraits<ExtField>::setPrimitive(getPrimitive<ExtField>());

    const int Dim = 2;
    const int r = 2; // field F_q, where q = r^2 -- we have F_4, so r = 2
    const int n = 8;
    typedef BMSDecoding<Dim, HermitianCodeParams<r, ExtField> > BMSDecoderT;
    BMSDecoderT bms_decoder(5); // C_4 code
    BMSDecoderT::FieldElemsCollection e;
    e.resize(n);

    // error positions!
    int pos1 = 1, pos2 = 7;
    e[pos1] = FieldElemTraits<ExtField>::multId(); // = 1
    e[pos2] = e[pos1];

    auto locs = bms_decoder.decode(e);
    auto refLocs = decltype(locs){pos1, pos2};
    ASSERT_EQUAL(locs, refLocs);
}

void runSuites() {
    cute::ide_listener</* empty for no IDE listener in standalone CUTE 2 */> lis;

    cute::suite PolyIOSuite;
    PolyIOSuite.push_back(CUTE(ioTestFor1Poly));
    PolyIOSuite.push_back(CUTE(outputTest));
    PolyIOSuite.push_back(CUTE(inputTestForNPoly));
    PolyIOSuite.push_back(CUTE(inputTestForNPolyOverGF));
    PolyIOSuite.push_back(CUTE(polySubscript));
    PolyIOSuite.push_back(CUTE(testPolyToDegCoefMapConversion));
    PolyIOSuite.push_back(CUTE(polyPowerPrinting));
    cute::makeRunner(lis)(PolyIOSuite, "The Polynomial Input-Output Suite");

    cute::suite PointSuite;
    PointSuite.push_back(CUTE(pointComparison));
    PointSuite.push_back(CUTE(pointIncreasing));
    PointSuite.push_back(CUTE(pointCollectionOperations));

    cute::suite PolynomialArithmeticSuite;
    PolynomialArithmeticSuite.push_back(CUTE(convolutionTest));
    PolynomialArithmeticSuite.push_back(CUTE(monomialMultiplication));
    PolynomialArithmeticSuite.push_back(CUTE(scalarMultiplication));
    PolynomialArithmeticSuite.push_back(CUTE(summation));
    PolynomialArithmeticSuite.push_back(CUTE(equality));
    PolynomialArithmeticSuite.push_back(CUTE(eval));

    cute::suite bmsaTestingSuite;
    bmsaTestingSuite.push_back(CUTE(sakatasExample2D));
    bmsaTestingSuite.push_back(CUTE(sakatasExample3D));

    cute::suite bmsaDecoding;
    bmsaDecoding.push_back(CUTE(curveArithmetic));
    bmsaDecoding.push_back(CUTE(bmsaDecodingCLOS05Example));

    cute::makeRunner(lis)(PointSuite, 
            "The Point Suite");
    cute::makeRunner(lis)(PolynomialArithmeticSuite,
            "The Polynomial Arithmetic Suite");
    cute::makeRunner(lis)(bmsaTestingSuite,
            "The BMS-algorithm Testing Suite");
    cute::makeRunner(lis)(bmsaDecoding,
            "The BMS-algorithm-based Decoding Algorithm Testing Suite");
}

}  // namespace TestMVPoly

int main() {
    google::InitGoogleLogging("mv-poly");
    google::InstallFailureSignalHandler();
    LOG(INFO) << "Let the tests start!\n";
    TestMVPoly::runSuites();
}
