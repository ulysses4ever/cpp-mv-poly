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
    /// Creates point (0, 0, ..., 0).
    Point() {
        data.assign(0);
    }

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
    bool operator==(Point<Dim> const & lhs, Point<Dim> const & rhs) {
        return lhs.data == rhs.data;
    }

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
 * @return Result of comparison two points by current monomial order. True is
 * \c lhs less then \c rhs, false otherwise.
 */
template<int Dim>
inline
bool totalLess(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return lhs.weight() < rhs.weight()
            || std::lexicographical_compare(lhs.rbegin(), lhs.rend(),
                    rhs.rbegin(), rhs.rend());
}

template<int Dim>
inline
bool operator<(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return totalLess(lhs, rhs);
}

/**
 * Reflexive version of \c totalLess.
 * @param[in] lhs Left-hand side argument of “less-or-equal”.
 * @param[in] rhs Right-hand side argument of “less-or-equal”.
 * @return Result of comparison two points by current monomial order. True is
 * \c lhs less then or equal to \c rhs, false otherwise.
 */
template<int Dim>
inline
bool totalLessOrEqual(Point<Dim> const & lhs, Point<Dim> const & rhs) {
    return totalLess(lhs, rhs) || (lhs == rhs);
}

/**
 * Point's pre-increment following antilexicographic order.
 */
template<int Dim>
inline
Point<Dim>& Point<Dim>::operator++() {
    typename Point<Dim>::reverse_iterator rit = this->rbegin();
    typename Point<Dim>::reverse_iterator rit_to_inc(rit++);
    while (rit != this->rend()) {
        // loop invariant: rit is one step further then rit_to_inc
        if (*rit > 0) {
            --*rit;
            ++*rit_to_inc;
            return *this;
        }
        ++rit;
        ++rit_to_inc;
    }
    // rit == pt.rend()
    // pt = (0, 0, ..., 0, a) and result shoud be (a + 1, 0, ..., 0)
    *rit_to_inc = *(this->rbegin()) + 1;
    *(this->rbegin()) = 0;
    return *this;
}

/**
 * Point's post-increment following antilexicographic order.
 */
template<int Dim>
inline
Point<Dim> Point<Dim>::operator++(int) {
    Point<Dim> old(*this);
    ++*this;
    return old;
}

/**
 * Coordinate-wise point summation.
 */
template<int Dim>
inline
Point<Dim> operator+(Point<Dim> lhs, Point<Dim> const & rhs) {
    // lhs: pass-by-copy optimization
    return lhs += rhs;
}

/**
 * Coordinate-wise point subtaraction.
 */
template<int Dim>
inline
Point<Dim> operator-(Point<Dim> lhs, Point<Dim> const & rhs) {
    return lhs -= rhs;
}

/**
 * Point slice. It is kind of Decorator (cf. [GoF]) for point instance which
 * shifts the index used in \c Point subscript operator on \c Offset
 * positions. It is used when subscript operation goes deeper in the nestedness
 * level of or main recursive type \c MVPolyType (which is
 * <tt>Polynomial<… Polynomial<T>… ></tt>).
 *
 * We probaly should use some polymorphism (either static or dynamic)
 * to allow \c pt field contain either Point or Slice itself, but we decided to
 * “optimise away” such possibility always keeping a Point instance — this is
 * managed by \c make_slice utility which we consider a part of Slice implementation
 * so this abstraction is not that leaky. As a proof of this we recall generic
 * utility \c apply_subscript which polymorphically calls operator[] for
 * either Point or Slice instance.
 *
 * @param Dim Initial dimension of the point beeing sliced.
 * @param Offset Starting index in initial point to subscript from in
 * the current slice.
 */
template<int Dim, int Offset>
class Slice {
    Point<Dim> const & pt;

public:
    Slice(Point<Dim> const & pt_) : pt(pt_) {}

    Point<Dim> const & getImpl() {return pt;}

    typedef typename Point<Dim>::const_reference const_reference;

    typedef typename Point<Dim>::size_type size_type;

    const_reference
    operator[](size_type n) const {
        return pt[n + Offset];
    }
};

/// Slice of the point is “1-slice”. For 1-slice sl[i] ~ pt[i + 1].
template<int Dim>
Slice<Dim, 1> make_slice(Point<Dim> const & pt) {
    return Slice<Dim, 1>(pt);
}

/** Slice of the Slice<Dim, Offset> is Slice<Dim, Offset + 1>. It is convinient
 * for us to slice by 1 position each time.
 */
template<int Dim, int Offset>
Slice<Dim, Offset + 1> make_slice(Slice<Dim, Offset> const & sl) {
    return Slice<Dim, Offset + 1>(sl.getImpl());
}

#endif /* POINT_HPP_ */
