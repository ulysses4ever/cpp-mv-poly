/*
 * Point.hpp
 *
 *  Created on: 13.07.2010
 *      Author: ulysses
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

};

template<int Dim>
inline
std::ostream& operator<<(std::ostream & os, Point<Dim> const & p) {
    std::copy(p.begin(), p.end(), std::ostream_iterator<int>(os, " "));
    return os;
}

/**
 * Compare two points by coordinates. <tt>byCoordinateLess(lhs, rhs) == true</tt> iff
 * lhs[i] <= rhs[i] for all i. It is a partial order on the set of points.
 * @param lhs Left-hand side argument of “less”.
 * @param rhs Right-hand side argument of “less”.
 * @return Result of comparison two points by coordinates: true if
 * lhs[i] <= rhs[i] for all i false otherwise.
 */
template<int Dim>
inline
bool byCoordinateLess(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), true,
            std::logical_and<bool>(), std::less_equal<int>());
}


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
