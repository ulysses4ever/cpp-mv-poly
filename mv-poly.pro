HEADERS += \
    Utilities.hpp \
    Point.hpp \
    NtlUtilities.hpp \
    NtlPolynomials.hpp \
    mv_poly.hpp \
    CurveArithmetic.hpp \
    CoefficientTraits.hpp \
    bmsa-decoding.hpp \
    bmsa.hpp \
    BasicExamples.hpp

SOURCES += \
    main.cpp

unix|win32: LIBS += -lntl

unix|win32: LIBS += -lglpk

unix|win32: LIBS += -lglog

QMAKE_CXXFLAGS += -std=c++0x
