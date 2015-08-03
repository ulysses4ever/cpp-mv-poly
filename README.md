# cpp-mv-poly

(Automatically exported from code.google.com/p/cpp-mv-poly)

In the project we create an implementation of multivariate polynomials arithmetic which will use modern template 
techniques, will be able to work with different fields (rings) for polynomial coefficients (e.g. we try to 
communicate with NTL (cf. [NTL]) data structures which represent those notions). As a part of our “technical 
research” on using template techniques we are interested in applying template programming design patterns (e.g. 
from [VJ]) as well as adapting well-known OOP design patterns descussed in [GoF] for the use in generic 
programming setting. 

The main goal of the project though not polynomials implementation itself but construction of 
BMS-algorithm (cf. [CLO'S, ch. 10]).

### References 

  * [NTL] _NTL: A Library for doing Number Theory_ by Victor Shoup, http://shoup.net/ntl/
  * [VJ] _D. Vandevoorde, N. M. Josuttis_, C++ Templates: The Complete Guide.
  * [GoF] _E. Gamma, R. Helm, R. Johnson, and J. M. Vlissides_, Design Patterns: Elements of Reusable Object-Oriented Software.
  * [CLO'S] _D. A. Cox, J. Little, D. O'Shea_, Using Algebraic Geometry.
