/*
 * utilities.hpp
 *
 *  Created on: 21.06.2011
 *      Author: ulysses
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <sstream>
#include <string>
#include <utility>
#include <iostream>

#define ARR_LEN(arr) (sizeof(arr) / sizeof(*(arr)))

//# if defined ERROR
//#   undef ERROR
//# endif
//# define ERROR() { fprintf(stderr, "Error\n"); exit(1); }

namespace mv_poly {

template<typename T1, typename T2>
inline
std::ostream& operator<<(std::ostream & os, std::pair<T1, T2> const & p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

template<typename M>
std::string mapToStr(M const & m) {
    std::ostringstream oss;
    for (auto pr : m) 
        oss
              << pr.first  << " : "
              << pr.second << " ";
}

} // namespace mv_poly

#endif /* UTILITIES_HPP_ */
