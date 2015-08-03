# cpp-mv-poly

(Automatically exported from code.google.com/p/cpp-mv-poly)

In the project we create an implementation of multivariate polynomials arithmetic which use modern template 
techniques, is be able to work with different fields (rings) for polynomial coefficients (e.g. we used to apply NTL (cf. [NTL]) data structures which represent those notions). As a part of our “technical 
research” on using template techniques we are interested in applying template programming design patterns (e.g. 
from [VJ]) as well as adapting well-known OOP design patterns descussed in [GoF] for the objects from abstract algebra. More precisely the main goal of the project is construction of BMS-algorithm (cf. [CLO'S, ch. 10]).

### Dependencies

  * [NTL 5+](http://shoup.net/ntl/): A Library for doing Number Theory;
  * [CUTE 2+](http://cute-test.com/projects/cute/wiki/CUTE_standalone): “C++ Unit Testing Easier”;
  * [GLPK](http://www.gnu.org/software/glpk/): GNU Linear Programming Kit;
  * [google-glog](http://code.google.com/p/google-glog/): Logging library for C++ (when logging is turned off, library is not needed for building a project).

### Build

The project is headers-only, so in order to try it out you have to use some test source code (in .cpp-file). Such a role is played by Test.cpp from a repo. So current version of project could be build e.g. via:

    g++ -std=c++11 -I/path/to/cute/cute_lib -o Test Test.cpp -lntl -lglpk

You need to have a compiler with thorough support for C++11 (say, GCC 4.7+). And the libraries listed in the _Dependencies_ section.

### References 

  * [NTL] _NTL: A Library for doing Number Theory_ by Victor Shoup, http://shoup.net/ntl/
  * [VJ] _D. Vandevoorde, N. M. Josuttis_, C++ Templates: The Complete Guide.
  * [GoF] _E. Gamma, R. Helm, R. Johnson, and J. M. Vlissides_, Design Patterns: Elements of Reusable Object-Oriented Software.
  * [CLO'S] _D. A. Cox, J. Little, D. O'Shea_, Using Algebraic Geometry.
