/**
 * Describes Point class template with accompanying routines.
 *
 * @date 2010-07-13
 * @author ulysses
 */

#ifndef POINT_HPP_
#define POINT_HPP_

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

    typename ImplType::reference
    operator[](typename ImplType::size_type n) {
        return data[n];
    }

    typename ImplType::const_reference
    operator[](typename ImplType::size_type n) const {
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
 * lhs[i] <= rhs[i] for all i false otherwise.
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
//    int lhsPow = lhs.weight();
//    int rhsPow = rhs.weight();
//    if (lhsPow == rhsPow)
//        return lhs[0] > rhs[0];
//    else
//        return lhsPow < rhsPow;
    return lhs.weight() < rhs.weight()
            || std::lexicographical_compare(lhs.rbegin(), lhs.rend(),
                    rhs.rbegin(), rhs.rend());
}

#endif /* POINT_HPP_ */
