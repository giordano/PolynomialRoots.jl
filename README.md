# CmplxRoots

[![Build Status](https://travis-ci.org/giordano/CmplxRoots.jl.svg?branch=master)](https://travis-ci.org/giordano/CmplxRoots.jl) [![Coverage Status](https://coveralls.io/repos/github/giordano/CmplxRoots.jl/badge.svg?branch=master)](https://coveralls.io/github/giordano/CmplxRoots.jl?branch=master) [![codecov.io](https://codecov.io/github/giordano/CmplxRoots.jl/coverage.svg?branch=master)](https://codecov.io/github/giordano/CmplxRoots.jl?branch=master) <!-- [![CmplxRoots](http://pkg.julialang.org/badges/CmplxRoots_0.4.svg)](http://pkg.julialang.org/?pkg=CmplxRoots) [![CmplxRoots](http://pkg.julialang.org/badges/CmplxRoots_0.5.svg)](http://pkg.julialang.org/?pkg=CmplxRoots) -->

Introduction
------------

`CmplxRoots.jl` is a library for finding roots of complex polynomials, written
in [Julia](http://julialang.org/).

This is just a Julia wrapper around the
[General Complex Polynomial Root Solver](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/),
written in Fortran, by **Jan Skowron** and **Andy Gould**.  All the credits goes
to them for the underlying functions, blame me for any problem with the Julia
interface.  Feel free to report bugs and make suggestions at
https://github.com/giordano/CmplxRoots.jl/issues.

The root finding algorithm employed in this library is described in

* J. Skowron & A. Gould, 2012, "General Complex Polynomial Root Solver and Its
  Further Optimization for Binary Microlenses",
  [arXiv:1203.1034](http://arxiv.org/abs/1203.1034)

This algorithm aims to be fast and precise, more than the well known `zroots`
procedure described in *Numerical Recipes* books, whose implementations in C and
Fortran are not available as free software, according to the
[definition](https://www.gnu.org/philosophy/free-sw.html) of the Free Software
Foundation.

**Note:** This package has been tested only on GNU/Linux and OS X systems.
Trying to install it on Windows will likely fail, please report at
https://github.com/giordano/CmplxRoots.jl/issues/1 if you manage to install the
package on this system.

Installation
------------

`CmplxRoots.jl` is available for Julia 0.4 and later versions, and can be installed
with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.clone("https://github.com/giordano/CmplxRoots.jl.git")
```

Installation script will download General Complex Polynomial Root Solver source
code and build the shared object.  In order to accomplish this task a Fortran
compiler is needed.

<!-- You may need to update your package list with `Pkg.update()` in order to get the -->
<!-- latest version of `CmplxRoots.jl`. -->

Usage
-----

After installing the package, run

``` julia
using CmplxRoots
```

or put this command into your Julia script.

`CmplxRoots.jl` provides two functions to find the roots of complex polynomials

``` julia
roots(polynomial[, roots, polish=true])
roots5(polynomial[, roots, polish=true])
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
  be complex or real

Optional arguments are:

* `roots`: if you roughly know in advance the position of the roots, you can
  pass the vector with the known roots to speed up convergence.  `roots` vector
  must be one element shorther than `polynomial`.  Elements of the vector can be
  complex or real
* `polish` (boolean keyword)

	* in `roots`: if set to `true`, after all roots have been found by dividing
	  original polynomial by each root found, all roots will be polished using
	  full polynomial.  Default is `false`
	* in `roots5`: if you know the roots pretty well, for example because you
	  have changed the coefficiens of the polynomial only a bit, so the two
	  closest roots are most likely still the closest ones, set this keyword to
	  `true`.  Default it `false`

The functions return in output the `Complex128` vector with all roots of
`polynomial`.  **Note:** if `roots` optional argument is provided, it is *not*
changed in-place.

Example
-------

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
         -0.0+5.0im
 4.0-1.11022e-15im
  1.51244e-16+3.0im
  2.0+2.6822e-16im
 -8.97111e-18+1.0im

julia> roots5([120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1])
5-element Array{Complex{Float64},1}:
          0.0+5.0im
 4.0-5.92119e-16im
 2.0-9.08596e-17im
  2.66454e-15+3.0im
 -2.66454e-15+1.0im
```

`CmplxRoots.jl` handles polynomials with high-multiplicity roots as well.  For
example, consider

```
(x - 1)^5 = x^5 - 5x^4 + 10x^3 - 10x^2 + 5x -1
```

You can find its roots with

``` julia
julia> roots([1, -5, 10, -10, 5, -1])
5-element Array{Complex{Float64},1}:
 1.0-0.0im
 1.0+0.0im
 1.0+0.0im
 1.0+0.0im
 1.0+0.0im

julia> roots5([1, -5, 10, -10, 5, -1])
5-element Array{Complex{Float64},1}:
 1.0+0.0im
 1.0+0.0im
 1.0+0.0im
 1.0+0.0im
 1.0-0.0im
```

Related projects
----------------

Another Julia package for finding roots of complex polynomials is
[`Polynomials.jl`](https://github.com/Keno/Polynomials.jl), by Jameson Nash and
other contributors.  Differently from `CmplxRoots.jl`, this is entirely written
in Julia and does not depend on packages written in foreign languages, so it
runs smoothly on all operating systems where Julia is available.  In order to
solve the polynomial, `Polynomials.jl` calculates eigenvalues of its companion
matrix.  This package does much more than finding roots of polynomials (among
other features, it can integrate and differentiate polynomials), but
`CmplxRoots.jl` is usually faster than `Polynomials.jl` and often slightly more
precise.  If you are after speed and precision, `CmplxRoots.jl` can still be
your choice.

License
-------

The `CmplxRoots.jl` package is licensed under the GNU Lesser General Public
License version 3 or any later version, as well as under a "customary scientific
license", which implies that if this code was important in the scientific
process or for the results of your scientific work, you are asked for the
appropriate citation of the paper Skowron & Gould 2012
(http://arxiv.org/abs/1203.1034).  The original author of `CmplxRoots.jl` is
Mosè Giordano.
