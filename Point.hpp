/** @file
 *
 * Describes Point class template with accompanying routines.
 *
 * @date 2010-07-13
 * @author ulysses
 */

#ifndef POINT_HPP_
#define POINT_HPP_

#include <tr1/array>

/**
 * Point in N-dimensional integer lattice.
 * @param Dim Dimension of point lattice.
 */
template<int Dim>
class Point {
    typedef std::tr1::array<int, Dim> ImplType;

    ImplType data;

public:
    int weight() const {
        return std::accumulate(data.begin(), data.end(), 0);
    }

    typedef typename ImplType::size_type size_type;

    typedef typename ImplType::reference reference;

    typedef typename ImplType::const_reference const_reference;

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
};

/**
 * Simple Point output.
 * @param[out] os Target output stream.
 * @param[in] pt Point to be outputed in \c os.
 * @return Output stream \c os after point \c p have been
 * outputed to \c os (conventionally).
 */
template<int Dim>
inline
std::ostream& operator<<(std::ostream & os, Point<Dim> const & pt) {
    std::copy(pt.begin(), pt.end(), std::ostream_iterator<int>(os, " "));
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
template<int Dim>
inline
bool byCoordinateLess(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), true,
            std::logical_and<bool>(), std::less_equal<int>());
}

/**
 * This total oreder predicate implements monomial order. Namely graded
 * antilexicographic order.
 * @param[in] lhs Left-hand side argument of “less”.
 * @param[in] rhs Right-hand side argument of “less”.
 * @return Result of comparison two points by current monomial order.
 */
template<int Dim>
inline
bool totalLess(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return lhs.weight() < rhs.weight()
            || std::lexicographical_compare(lhs.rbegin(), lhs.rend(),
                    rhs.rbegin(), rhs.rend());
}

/**
 * In place point increasing following antilexicographic order.
 * @param[in,out] pt Point to be increased.
 */
template<int Dim>
inline
void increase(Point<Dim> & pt) {
    typename Point<Dim>::reverse_iterator rit = pt.rbegin();
    typename Point<Dim>::reverse_iterator rit_to_inc(rit++);
    while (rit != pt.rend()) {
        // loop invariant: rit is one step further then rit_to_inc
        if (*rit > 0) {
            --*rit;
            ++*rit_to_inc;
            return;
        }
        ++rit;
        ++rit_to_inc;
    }
    // rit == pt.rend()
    // pt = (0, 0, ..., 0, a) and result shoud be (a + 1, 0, ..., 0)
    *rit_to_inc = *pt.rbegin() + 1;
    *pt.rbegin() = 0;
}

/**
 * Point slice.
 */
template<typename Pt>
class Slice {
    Pt const & pt;

public:
    Slice(Pt const & pt_) : pt(pt_) {}

//    typename ImplType::reference
//    operator[](typename ImplType::size_type n) {
//        return pt[n + Beg];
//    }
    typedef typename Pt::const_reference const_reference;

    typedef typename Pt::size_type size_type;

    const_reference
    operator[](size_type n) const {
        return pt[n + 1];
    }
};

template<int Dim>
Point<Dim - 1> make_slice(Point<Dim> const & pt) {
    Point<Dim - 1> result;
    for (int i = 0; i < Dim - 1; ++i)
        result[i] = pt[i+1];
    return result;
}

#endif /* POINT_HPP_ */
