/** @file Point.hpp
 *
 * Describes Point class template with accompanying routines.
 *
 * @date 2010-07-13
 * @author ulysses
 */

#ifndef POINT_HPP_
#define POINT_HPP_

#include <array>
#include <algorithm>
#include <iterator>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <cstdio>

#include <tr1/functional>

#include <boost/foreach.hpp>
//#include <boost/lexical_cast.hpp>

//#include <lpsolve/lp_lib.h>
#include <glpk.h>

#include "Utilities.hpp"

namespace mv_poly {

template<typename PointImpl>
struct GradedAntilexMonomialOrder;

template<
    int Dim,
    template <typename> class OrderPolicy
        = GradedAntilexMonomialOrder
>
class Point;

template<int Dim, template <typename> class OrderPolicy>
bool operator==(
        Point<Dim, OrderPolicy> const & lhs,
        Point<Dim, OrderPolicy> const & rhs);

/**
 * \class Point
 * Point in N-dimensional integer lattice.
 * @param Dim Dimension of point lattice.
 */
template<
    int Dim,
    template <typename PointImpl> class OrderPolicy
//        = GradedAntilexMonomialOrder
> class Point : OrderPolicy< std::array<long, Dim> > {

    typedef typename Point::PointImplType ImplType; // PointImplType inherited
                                            // from GradedAntilexMonomialOrder
                                            // to avoid duplication

    ImplType data;

public:
    /// Creates point (0, 0, ..., 0).
    Point() {
        data.fill(0);
    }

    bool operator<(Point<Dim, OrderPolicy> const & other) const {
        return totalLess(data, other.data);
    }

    typedef typename ImplType::value_type value_type;

    typedef typename ImplType::size_type size_type;

    typedef typename ImplType::reference reference;

    typedef typename ImplType::const_reference const_reference;

    Point(std::initializer_list<value_type> data) {
        std::copy(data.begin(), data.end(), this->data.begin());
    }

    reference
    operator[](size_type n) {
        return data[n];
    }

    const_reference
    operator[](size_type n) const {
        return data[n];
    }

    typedef typename ImplType::iterator iterator;

    typedef typename ImplType::const_iterator const_iterator;

    iterator
    begin()
    { return data.begin(); }

    const_iterator
    begin() const
    { return data.begin(); }

    iterator
    end()
    { return data.end(); }

    const_iterator
    end() const
    { return data.end(); }

    typedef typename ImplType::reverse_iterator reverse_iterator;

    typedef typename ImplType::const_reverse_iterator const_reverse_iterator;

    reverse_iterator
    rbegin()
    { return data.rbegin(); }

    const_reverse_iterator
    rbegin() const
    { return data.rbegin(); }

    reverse_iterator
    rend()
    { return data.rend(); }

    const_reverse_iterator
    rend() const
    { return data.rend(); }

    Point& operator+=(Point const & other) {
        std::transform(this->begin(), this->end(),
                other.begin(), this->begin(), std::plus<int>());
        return *this;
    }

    Point& operator-=(Point const & other) {
        std::transform(this->begin(), this->end(),
                other.begin(), this->begin(), std::minus<int>());
        return *this;
    }

    Point& operator++();

    Point operator++(int);

    /**
     * Compaison for equality: two points are equal iff all corresponding coordinates
     * are equal.
     * @param lhs Left-hand side argument for testing on equality.
     * @param rhs Right-hand side argument for testing on equality.
     * @return True if all corresponding coordinates of \c lhs and \c rhs are equal,
     * false otherwise.
     */
    friend
    bool operator==<>(
            Point<Dim, OrderPolicy> const & lhs,
            Point<Dim, OrderPolicy> const & rhs);
;

}; // class Point

template<int Dim, template <typename> class OrderPolicy>
bool operator==(
        Point<Dim, OrderPolicy> const & lhs,
        Point<Dim, OrderPolicy> const & rhs) {
   return lhs.data == rhs.data;
}

/**
 * Summs up elements of a given container \c c.
 * @param c Container to summ up elements.
 * @return Summa of the elements of \c c.
 */
template<typename Cont>
inline
int weight(Cont const & c) {
    return std::accumulate(c.begin(), c.end(), 0);
}

template<typename PointImpl>
struct GradedAntilexMonomialOrder {

    typedef PointImpl PointImplType;

    /**
     * This total oreder predicate implements graded antilexicographic order.
     * @param[in] lhs Left-hand side argument of “less”.
     * @param[in] rhs Right-hand side argument of “less”.
     * @return Result of comparison two points by current monomial order. True is
     * \c lhs less then \c rhs, false otherwise (rhs is less or equal to lhs).
     */
    static bool totalLess(PointImplType const & lhs, PointImplType const & rhs) {
        int lw = weight(lhs);
        int rw = weight(rhs);
        return (lw < rw)
                || (lw == rw
                        && std::lexicographical_compare(lhs.rbegin(), lhs.rend(),
                                rhs.rbegin(), rhs.rend()));
    }

    void inc(PointImplType & data) {
        using namespace std::tr1::placeholders;
        using std::tr1::bind;
        typename PointImplType::iterator
                itInc = std::find_if(data.begin(), data.end(),
                    std::tr1::bind(std::logical_not<bool>(),
                            std::tr1::bind(std::equal_to<int>(), 0, _1))),
                it(itInc++);
        if (it == data.end())
            data[0] = 1;
        else if (itInc == data.end()) {
            int a = *it + 1;
            *it = 0;
            data[0] = a;
        } else {
            ++(*itInc);
            int a = *it - 1;
            *it = 0;
            data[0] = a;
        }
    }
};

template <int a, int b>
class LpSolveHolder {
    static const int dim = 2;

    glp_prob *lp;

    glp_iocp parm;

public:
    LpSolveHolder() {
        glp_term_out(GLP_OFF);

        const int ne = 2;
        int ia[1 + ne], ja[1 + ne];
        double ar[1 + ne];

        lp = glp_create_prob();

        glp_set_obj_dir(lp, GLP_MIN);

        glp_add_rows(lp, 1);

        glp_add_cols(lp, 2);
        glp_set_col_bnds(lp, 1, GLP_DB, 0.0, b - 1); // x1
        glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0); // x2

        glp_set_col_kind(lp, 1, GLP_IV);
        glp_set_col_kind(lp, 2, GLP_IV);

        glp_set_obj_coef(lp, 1, a);
        glp_set_obj_coef(lp, 2, b);

        ia[1] = 1, ja[1] = 1, ar[1] = a;
        ia[2] = 1, ja[2] = 2, ar[2] = b;

        glp_load_matrix(lp, 2, ia, ja, ar);
        //glp_simplex(lp, NULL);
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
    }

    template <typename Cont>
    void lpInc(Cont & c) {
        glp_set_row_bnds(lp, 1, GLP_LO, c[0]*a + c[1]*b + 1, 0.0);
        glp_intopt(lp, &parm);
        c[0] = glp_mip_col_val(lp, 1);
        c[1] = glp_mip_col_val(lp, 2);
    }

    ~LpSolveHolder() {
        glp_delete_prob(lp);
    }
};

template<int a, int b>
struct WeightedOrder {

    template<typename PointImpl>
    struct impl {
        typedef PointImpl PointImplType;

        static bool totalLess(PointImplType const & lhs, PointImplType const & rhs) {
            return weight(lhs) < weight(rhs);
        }

        static void inc(PointImplType & data) {
            lp.lpInc(data);
        }


    private:
        static LpSolveHolder<a, b> lp;

        static int weight(PointImplType const & data) {
            return a*data[0] + b*data[1];
        }

    };

};

template <int a, int b>
template <typename T>
LpSolveHolder<a, b>
WeightedOrder<a, b>::impl<T>::lp = LpSolveHolder<a, b>();


/**
 * Simple Point output.
 * @param[out] os Target output stream.
 * @param[in] pt Point to be outputed in \c os.
 * @return Output stream \c os after point \c p have been
 * outputed to \c os (conventionally).
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
std::ostream& operator<<(std::ostream & os, Point<Dim, OrderPolicy> const & pt) {
    os << "(";
    //std::copy(pt.begin(), pt.end(), std::ostream_iterator<int>(os, ","));
    for(int i = 0; i < Dim - 1; ++i)
        os << pt[i] << ", ";
    os << pt[Dim - 1] << ")";
    return os;
}

/**
 * Simple 1-dimensional point output.
 * @param[out] os Target output stream.
 * @param[in] pt Point to be outputed in \c os.
 * @return Output stream \c os after point \c p have been
 * outputed to \c os (conventionally).
 */
template<template <typename PointImpl> class OrderPolicy>
inline
std::ostream& operator<<(std::ostream & os, Point<1, OrderPolicy> const & pt) {
    os << pt[0];
    return os;
}

/**
 * Compare two points by coordinates. <tt>byCoordinateLess(lhs, rhs) == true</tt> iff
 * lhs[i] <= rhs[i] for all i. It is a partial order on the set of points.
 * @param[in] lhs Left-hand side argument of “less”.
 * @param[in] rhs Right-hand side argument of “less”.
 * @return Result of comparison two points by coordinates: true if
 * lhs[i] <= rhs[i] for all i, false otherwise.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
bool byCoordinateLess(
        Point<Dim, OrderPolicy> const & lhs,
        Point<Dim, OrderPolicy> const & rhs) {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), true,
            std::logical_and<bool>(), std::less_equal<int>());
}

/**
 * Checks if given point \pt  is by-coordinate less then any point in
 * range given by iterators \c beg and \c end.
 * @param[in] pt Point to check for by-coordinate less.
 * @param[in] beg Iterator that points to the begining (inclusively) of the
 * point sequence to be check against.
 * @param[in] end Iterator that points to the end (exclusively) of the
 * point sequence to be check against.
 * @return True if there exists a point pt' in [*beg, *end) such that
 * byCoordinateLess(pt, pt'), false otherwise.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy, typename It>
inline
bool byCoordinateLessThenAny(Point<Dim, OrderPolicy> const & pt, It beg, It end) {
    using std::tr1::bind;
    using std::tr1::placeholders::_1;
    using std::tr1::cref;
    return end != std::find_if(beg, end,
            std::tr1::bind(&byCoordinateLess<Dim, OrderPolicy>, cref(pt), _1));
}

/**
 * Checks if given point \pt  is by-coordinate less then any point in
 * collection \c c.
 * @param[in] pt Point to check for by-coordinate less.
 * @param[in] c Container of points to be check against.
 * @return True if there exists a point pt' in \c c such that
 * byCoordinateLess(pt, pt'), false otherwise.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy, typename PtCont>
inline
bool byCoordinateLessThenAny(Point<Dim, OrderPolicy> const & pt, PtCont const & c) {
    return byCoordinateLessThenAny(pt, c.begin(), c.end());
}

/**
 * Checks if given point \pt  is by-coordinate greater then any point in
 * collection \c c.
 * @param[in] pt Point to check for by-coordinate greater.
 * @param[in] c Container of points to be check against.
 * @return True is there exists a point pt' in \c c such that
 * byCoordinateLess(pt, pt'), false otherwise.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy, typename PtCont>
inline
bool byCoordinateGreaterThenAny(Point<Dim, OrderPolicy> const & pt,
        PtCont const & c) {
    using std::tr1::bind;
    using std::tr1::placeholders::_1;
    using std::tr1::cref;
    return c.end() != std::find_if(c.begin(), c.end(),
            bind(&byCoordinateLess<Dim, OrderPolicy>, _1, cref(pt)));
}

/**
 * Reflexive version of \c operator<.
 * @param[in] lhs Left-hand side argument of “less-or-equal”.
 * @param[in] rhs Right-hand side argument of “less-or-equal”.
 * @return Result of comparison two points by current monomial order. True is
 * \c lhs less then or equal to \c rhs, false otherwise.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
bool operator<=(Point<Dim, OrderPolicy> const & lhs,
        Point<Dim, OrderPolicy> const & rhs) {
    return (lhs < rhs) || (lhs == rhs);
}

/**
 * Point's pre-increment.
 *
 * @note Use OrderPolicy.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
Point<Dim, OrderPolicy>&
Point<Dim, OrderPolicy>::operator++() {
    inc(data);
    return *this;
}

/**
 * Point's post-increment.
 *
 * @note Use OrderPolicy.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
Point<Dim, OrderPolicy> Point<Dim, OrderPolicy>::operator++(int) {
    Point<Dim, OrderPolicy> old(*this);
    ++*this;
    return old;
}

/**
 * Coordinate-wise point summation.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
Point<Dim, OrderPolicy>
operator+(Point<Dim, OrderPolicy> lhs, Point<Dim, OrderPolicy> const & rhs) {
    // lhs: pass-by-copy optimization
    return lhs += rhs;
}

/**
 * Coordinate-wise point subtaraction.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy>
inline
Point<Dim, OrderPolicy>
operator-(Point<Dim, OrderPolicy> lhs, Point<Dim, OrderPolicy> const & rhs) {
    return lhs -= rhs;
}

/**
 * Gets all partial maximums from collection of \c Point (\c points) with
 * respect to by-coordinate partial order (cf. \c byCoordinateLess).
 * @param points Collection of points to be looked through for the maximums.
 * @return Maximum points with respect to by-coordinate partial order
 * (cf. \c byCoordinateLess) from the \c points.
 */
template<
    int Dim,
    template <typename PointImpl> class OrderPolicy,
    template<typename T, typename S = std::allocator<T> > class Cont>
Cont<Point<Dim, OrderPolicy> >
getPartialMaximums(Cont<Point<Dim, OrderPolicy> > const & points) {
    using std::tr1::bind;
    using std::tr1::placeholders::_1; // usually we use “using” directive:
    // using namespace std::tr1::placeholders; — but here it yelds some ambiguity
    // while resolving _1 symbol, so we use declaration. Sad but true...
    using std::tr1::cref;
    using std::find_if;
    Cont<Point<Dim, OrderPolicy> > result;
    typedef Point<Dim, OrderPolicy> Pt;
    BOOST_FOREACH(Pt const & pt, points) {
        // if there is a point in result which dominates pt, then throw pt away
        if (byCoordinateLessThenAny(pt, result))
            continue;
        // else we add pt to the result but first delete all points
        // in result dominated by pt
        result.erase(
                remove_if(result.begin(), result.end(),
                        bind(&byCoordinateLess<Dim, OrderPolicy>, _1, cref(pt))),
                result.end());
        result.push_back(pt);
    }
    return result;
}

/**
 * Complimentary to \c getPartialMaximums.
 * @param points Collection of points to be looked through for the maximums.
 * @return Maximum points with respect to by-coordinate partial order
 * (cf. \c byCoordinateLess) from the \c points.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy,
    template<typename T, typename S = std::allocator<T> > class Cont>
Cont<Point<Dim, OrderPolicy> >
getPartialMinimums(Cont<Point<Dim, OrderPolicy> > const & points) {
    // TODO: probably we should refactor out this function together with
    // getPartialMaximums to get one generic function
    using std::tr1::bind;
    using std::tr1::placeholders::_1;
    using std::tr1::cref;
    using std::find_if;
    Cont<Point<Dim, OrderPolicy> > result;
    typedef Point<Dim, OrderPolicy> Pt;
    BOOST_FOREACH(Pt const & pt, points) {
        // if there is a point in result which dominates pt, then throw pt away
        if (byCoordinateGreaterThenAny(pt, result))
            continue;
        // else we add pt to the result but first delete all points
        // in result dominated by pt
        result.erase(
                remove_if(result.begin(), result.end(),
                        bind(&byCoordinateLess<Dim, OrderPolicy>, _1, cref(pt))),
                result.end());
        result.push_back(pt);
    }
    return result;
}

template<int Dim, template <typename PointImpl> class OrderPolicy,
    template<typename T, typename S = std::allocator<T> > class Cont>
Cont<Point<Dim, OrderPolicy> >
getConjugatePointCollection(Cont<Point<Dim, OrderPolicy> > const & points) {
    // construct some finite approximation set of Sigma-set, which contains conjugate
    // point set to be found; then we use exhaustive search in this finite set for
    // finding extremums (minimums in this case) as usual (see getPartialMaximums);
    // the approximation set defined as follows (and this definition is a
    // subject of probable future optimization, specification
    // and even correction ;)): it is all points totally
    // less (w.r.t. monomial order) then the point which has a weight equal to
    // the maximum weight of points in given collection plus 2
    // but the (totally) least from all such points
    // TODO: proof correctness of the algorithm
    using std::tr1::bind;
    using std::tr1::function;
    using std::tr1::placeholders::_1;
    using std::tr1::placeholders::_2;
    if (points.empty())
        return Cont< Point<Dim, OrderPolicy> >(1);
    //function<int (Point<Dim>)> pointWeight = bind(&Point<Dim>::weight, _1);
    int max_weight = weight(*std::max_element(points.begin(), points.end(),
            std::tr1::bind(std::less<int>(),
                    std::tr1::bind(&weight< Point<Dim, OrderPolicy> >, _1),
                    std::tr1::bind(&weight< Point<Dim, OrderPolicy> >, _2))));
    Point<Dim, OrderPolicy> upperPoint;
    upperPoint[0] = max_weight + 2;
    Cont< Point<Dim, OrderPolicy> > approxSigmaSet;
    for (Point<Dim, OrderPolicy> i; i < upperPoint; ++i) {
        if (! byCoordinateLessThenAny(i, points))
            approxSigmaSet.push_back(i);
    }
//    std::cout << "app set" << std::endl;
//    std::copy(approxSigmaSet.begin(), approxSigmaSet.end(),
//            std::ostream_iterator< Point<Dim> >(std::cout, "\n"));
    return getPartialMinimums(approxSigmaSet);
}
/**
 * Point slice. It is kind of Decorator (cf. [GoF]) for point instance which
 * shifts the index used in \c Point subscript operator on \c Offset
 * positions. It is used when subscript operation goes deeper in the nestedness
 * level of or main recursive type \c MVPolyType (which is
 * <tt>Polynomial<… Polynomial<T>… ></tt>).
 *
 * We probaly should use some polymorphism (either static or dynamic)
 * to allow \c pt field contain either Point or ConstSlice itself, but we decided to
 * “optimise away” such possibility always keeping a Point instance — this is
 * managed by \c make_slice utility which we consider a part of ConstSlice implementation
 * so this abstraction is not that leaky. As a proof of this we recall generic
 * utility \c apply_subscript which polymorphically calls operator[] for
 * either Point or ConstSlice instance.
 *
 * @param Dim Initial dimension of the point beeing sliced.
 * @param Offset Starting index in initial point to subscript from in
 * the current slice.
 */
template<typename Body, int Dim, int Offset>
class ConstSlice {
    Body const & pt;

public:
    ConstSlice(Body const & pt_) : pt(pt_) {}

    Body const & getImpl() const {return pt;}

    operator Body const &() {return pt;}

    typedef typename Body::const_reference const_reference;

    typedef typename Body::size_type size_type;

    const_reference
    operator[](size_type n) const {
        return pt[n + Offset];
    }

};

/// ConstSlice of the point is “1-slice”, for 1-slice sl[i] ~ pt[i + 1].
template<int Dim, typename Body>
ConstSlice<Body, Dim, 1>
make_slice(Body const & pt) {
    return ConstSlice<Body, Dim, 1>(pt);
}

/** ConstSlice of the ConstSlice<Dim, Offset> is ConstSlice<Dim, Offset + 1>.
 * It is convinient for us to slice by 1 position each time.
 */
template<typename Body, int Dim, int Offset>
ConstSlice<Body, Dim, Offset + 1>
make_slice(ConstSlice<Body, Dim, Offset> const & sl) {
    return ConstSlice<Body, Dim, Offset + 1>(sl.getImpl());
}

/**
 * Mutable version of ConstSlice.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy, int Offset>
class Slice {
    Point<Dim, OrderPolicy> & pt;

public:
    Slice(Point<Dim, OrderPolicy> & pt_) : pt(pt_) {}

    Point<Dim, OrderPolicy> & getImpl() const {return pt;}

    operator Point<Dim, OrderPolicy> &() {return pt;}

    typedef typename Point<Dim, OrderPolicy>::const_reference const_reference;

    typedef typename Point<Dim, OrderPolicy>::reference reference;

    typedef typename Point<Dim, OrderPolicy>::size_type size_type;

    const_reference
    operator[](size_type n) const {
        return pt[n + Offset];
    }

    reference
    operator[](size_type n) {
        return pt[n + Offset];
    }
};

/// ConstSlice of the point is “1-slice”, for 1-slice sl[i] ~ pt[i + 1].
template<int Dim, template <typename PointImpl> class OrderPolicy>
Slice<Dim, OrderPolicy, 1>
make_slice(Point<Dim, OrderPolicy> & pt) {
    return Slice<Dim, OrderPolicy, 1>(pt);
}

/** ConstSlice of the ConstSlice<Dim, Offset> is ConstSlice<Dim, Offset + 1>.
 * It is convinient for us to slice by 1 position each time.
 */
template<int Dim, template <typename PointImpl> class OrderPolicy, int Offset>
Slice<Dim, OrderPolicy, Offset + 1>
make_slice(Slice<Dim, OrderPolicy, Offset> const & sl) {
    return Slice<Dim, OrderPolicy, Offset + 1>(sl.getImpl());
}


} // namespace mv_poly
#endif /* POINT_HPP_ */
