# PolynomialRoots

[![Travis Build Status on GNU/Linux and OS X](https://travis-ci.org/giordano/PolynomialRoots.jl.svg?branch=master)](https://travis-ci.org/giordano/PolynomialRoots.jl) [![Appveyor Build Status on Windows](https://ci.appveyor.com/api/projects/status/jfa9e54lv92rqd3m?svg=true)](https://ci.appveyor.com/project/giordano/polynomialroots-jl) [![Coverage Status](https://coveralls.io/repos/github/giordano/PolynomialRoots.jl/badge.svg?branch=master)](https://coveralls.io/github/giordano/PolynomialRoots.jl?branch=master) [![codecov.io](https://codecov.io/github/giordano/PolynomialRoots.jl/coverage.svg?branch=master)](https://codecov.io/github/giordano/PolynomialRoots.jl?branch=master) [![PolynomialRoots](http://pkg.julialang.org/badges/PolynomialRoots_0.4.svg)](http://pkg.julialang.org/?pkg=PolynomialRoots) [![PolynomialRoots](http://pkg.julialang.org/badges/PolynomialRoots_0.5.svg)](http://pkg.julialang.org/?pkg=PolynomialRoots)

Introduction
------------

`PolynomialRoots.jl` is a library for finding roots of complex polynomials,
written in [Julia](http://julialang.org/).

This is an implementation in Julia of the
[General Complex Polynomial Root Solver](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/),
written in Fortran, by **Jan Skowron** and **Andrew Gould**.  All the credits
goes to them for the original algorithm.

The root finding algorithm employed in this library is described in

* J. Skowron & A. Gould, 2012, "General Complex Polynomial Root Solver and Its
  Further Optimization for Binary Microlenses",
  [arXiv:1203.1034](http://arxiv.org/abs/1203.1034)

This algorithm aims to be fast and precise, more than the well known `zroots`
procedure described in *Numerical Recipes* books, whose implementations in C and
Fortran are not available as free software, according to the
[definition](https://www.gnu.org/philosophy/free-sw.html) of the Free Software
Foundation.

`PolynomialRoots.jl` can also take advantage of native arbitrary precision
capabilities of Julia and achieve more precise results.

**Note**: the adopted algorithm can give highly inaccurate results for
polynomials of order above ~30.  This can be mitigated by using multiple
precision calculations (see example below).

Installation
------------

`PolynomialRoots.jl` is available for Julia 0.4 and later versions, and can be
installed with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.add("PolynomialRoots")
```

You may need to update your package list with `Pkg.update()` in order to get the
latest version of `PolynomialRoots.jl`.

Usage
-----

After installing the package, run

``` julia
using PolynomialRoots
```

or put this command into your Julia script.

`PolynomialRoots.jl` provides two functions to find the roots of complex
polynomials

``` julia
roots(polynomial[, roots, polish=true, epsilon=1e-20])
roots5(polynomial[, roots, epsilon=1e-20])
```

`roots` can be used to solve polynomials of any degree, `roots5` is tailored to
(and works only for) polynomials of degree 5.  This function exists because
[one of the methods](http://dx.doi.org/10.1086/309566) to calculate
[gravitational microlensing](https://en.wikipedia.org/wiki/Gravitational_microlensing)
amplification by a binary lens requires solving a fifth-order complex
polynomial.  `roots5` should be more robust than `roots` for this class of
polynomials.

The mandatory argument for both functions is:

* `polynomial`, the vector holding coefficients of the polynomial you want to
  solve, in ascending order, from the lowest to the highest.  Coefficients can
  be complex and real

Optional arguments are:

* `roots`: if you roughly know in advance the position of the roots, you can
  pass the vector with the known roots to speed up convergence.  `roots` vector
  must be one element shorther than `polynomial`.  In `roots5`, the roots will
  be only polished.  Elements of the vector can be complex and real
* `polish` (boolean keyword, only for `roots`): if set to `true`, after all
  roots have been found by dividing original polynomial by each root found, all
  roots will be polished using full polynomial.  Default is `false`
* `epsilon` (floating point optional keyword): this is used to determine the
  stopping criterion described in Skowron & Gould paper.  If not set, it
  defaults to machine precision of `polynomial` (and `roots`) argument(s).  This
  is *not* the precision with which the roots will be calculated.

The functions return in output the `Complex` vector with all roots of
`polynomial`.  **Note:** if `roots` optional argument is provided, it is *not*
changed in-place.

Examples
--------

Find the roots of

```
(x - i)*(x - 2)*(x - 3*i)*(x - 4)*(x - 5*i) =
  x^5 - (9*i + 6)*x^4 + (54*i - 15)*x^3 + (138 - 57*i)*x^2 - (184 + 90*i)*x + 120*i
```

This is a fifth-order polynomial, so we can find its zeros with both `roots` and
`roots5`:

``` julia
julia> roots([120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1])
5-element Array{Complex{Float64},1}:
 -1.55431e-15+5.0im
 4.0+1.77636e-16im
  1.55028e-15+3.0im
 -1.04196e-16+1.0im
 2.0-2.00595e-16im

julia> roots5([120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1])
5-element Array{Complex{Float64},1}:
  5.92119e-16+5.0im
  4.0-1.4803e-16im
 2.0+1.88202e-16im
 -1.88738e-15+3.0im
  2.10942e-15+1.0im
```

`PolynomialRoots.jl` handles polynomials with high-multiplicity roots as well.
For example, consider

```
(x + 1)^5 = x^5 + 5x^4 + 10x^3 + 10x^2 + 5x + 1
```

You can find its roots with

``` julia
julia> roots([1, 5, 10, 10, 5, 1])
5-element Array{Complex{Float64},1}:
 -1.0+0.0im
 -1.0+0.0im
 -1.0+0.0im
 -1.0+0.0im
 -1.0+0.0im

julia> roots5([1, 5, 10, 10, 5, 1])
5-element Array{Complex{Float64},1}:
 -1.0+0.0im
 -1.0+0.0im
 -1.0+0.0im
 -1.0+0.0im
 -1.0+0.0im
```

### Arbitrary precision ###

Due to limited precision of `Float64` type, extraction of roots of polynomials
can give inaccurate results, even for low-order polynomials.  This is caused,
i.e., by
[catastrophic cancellation](https://en.wikipedia.org/wiki/Loss_of_significance)
in computation of discriminant Δ = sqrt(b^2 - 4ac) of second-order polynomials.
[For example](http://www.cs.berkeley.edu/~wkahan/Qdrtcs.pdf), the actual roots
of

```
94906265.625*x^2 - 189812534*x + 94906268.375
```

are

```
x_1 = 1.000000028975958
x_2 = 1.000000000000000
```

but when you try to calculate them in double-precision you get

``` julia
julia> r = roots([94906268.375, -189812534, 94906265.625]);

julia> r[1]
1.0000000144879793 - 0.0im

julia> r[2]
1.0000000144879788 + 0.0im
```

If you are interested in double-precision accuracy, you can work around this
problem by calculating the roots with higher precision and then transforming the
result to double-precision.  Julia natively supports multiple precision
calculations, so what you have to do is only to pass `BigFloat` numbers to
`roots` function:

``` julia
julia> r = roots([BigFloat(94906268.375), BigFloat(-189812534), BigFloat(94906265.625)]);

julia> Float64(r[1])
1.0000000289759583

julia> Float64(r[2])
1.0
```

Note that in this case there is a trade-off between speed and higher accuracy
and precision.

Related projects
----------------

Another Julia package for finding roots of complex polynomials is
[`Polynomials.jl`](https://github.com/Keno/Polynomials.jl), by Jameson Nash and
other contributors.  This package does much more than finding roots of
polynomials (among other features, it can integrate and differentiate
polynomials).  In order to solve the polynomial, `Polynomials.jl` calculates
eigenvalues of its companion matrix, but `PolynomialRoots.jl` is usually faster
by up to an order of magnitude and often slightly more precise.  In addition,
`Polynomials` cannot extract roots in arbitrary precision.  If you are after
speed and precision, `PolynomialRoots.jl` can still be a better option.

How can I help?
---------------

Feel free to report bugs and sugges improvements at
https://github.com/giordano/PolynomialRoots.jl/issues.  You can also implement
other (possibly faster and/or more precise in more cases) root finding
algorithms and send a pull request at
https://github.com/giordano/PolynomialRoots.jl/pulls.

License
-------

The `PolynomialRoots.jl` package is licensed under version 2.0 of the Apache
License or the GNU Lesser General Public License version 3 or any later version,
at your option.  These are the same licenses used by the General Complex
Polynomial Root Solver.

A custom in the scientific comunity is (regardless of the licence you chose to
use or distribute this software under) that if this code was important in the
scientific process or for the results of your scientific work, we kindly ask you
for the **appropriate citation** of the paper Skowron & Gould 2012
(http://arxiv.org/abs/1203.1034), and we would be greatful if you pass the
information about the proper citation to anyone whom you redistribute this
software to.

The authors of the General Complex Polynomial Root Solver, the original Fortran
library (http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/) from which
`PolynomialRoots.jl` has been translated, are Jan Skowron, Andrew Gould.

The original author of `PolynomialRoots.jl` is Mosè Giordano.
