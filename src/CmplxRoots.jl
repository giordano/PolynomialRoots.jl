### CmplxRoots.jl --- Complex Polynomial Root Solver

# Copyright (C) 2016  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: polynomials, root finding

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

### Code:

module CmplxRoots

export roots, roots5

# Note: don't use Pkg.dir("PkgName") here because the package may be installed
# elsewhere.
const libcmplxroots = joinpath(dirname(@__FILE__), "..", "deps", "libcmplxroots")

function roots!(roots::Vector{Complex128}, poly::Vector{Complex128},
                degree::Integer, polish::Bool, use_roots::Bool)
    ccall((:cmplx_roots_gen_, libcmplxroots), Ptr{Void},
          (Ptr{Complex128}, # roots
           Ptr{Complex128}, # poly
           Ptr{Cint}, # degree
           Ptr{Cint}, # polish_roots_after
           Ptr{Cint}),# use_roots_as_starting_points
          roots, poly, Ref{Cint}(degree),
          Ref{Cint}(polish), Ref{Cint}(use_roots))
    return roots
end

"""
    roots(polynomial[, roots, polish=true]) -> roots

Find all the roots of `polynomial`, of any degree.

Arguments:

* `polynomial`: vector of coefficients (type `Number`) of the polynomial of
  which to find the roots, from the lowest coefficient to the highest one
* `roots` (optional argument): vector of initial guess roots.  If you have a
  very rough idea where some of the roots can be, this vector is used as
  starting value for Laguerre's method
* `polish` (optional boolean keyword): if set to `true`, after all roots have
  been found by dividing original polynomial by each root found, all roots will
  be polished using full polynomial.  Default is `false`

Function `root5` is specialized for polynomials of degree 5.
"""
function roots{N1<:Number,N2<:Number}(poly::Vector{N1}, roots::Vector{N2};
                                      polish::Bool=false)
    degree = length(poly) - 1
    @assert degree == length(roots) "`poly' must have one element more than `roots'"
    roots!(float(complex(roots)), float(complex(poly)), degree, polish, true)
end

function roots{N1<:Number}(poly::Vector{N1}; polish::Bool=false)
    degree = length(poly) - 1
    roots!(Array(Complex128, degree), float(complex(poly)), degree, polish, false)
end

function roots5!(roots::Vector{Complex128}, poly::Vector{Complex128},
                 polish::Bool)
    order_changed = Cint(0)
    ccall((:cmplx_roots_5_, libcmplxroots), Ptr{Void},
          (Ptr{Complex128}, # roots
           Ref{Cint}, # first_3_roots_order_changed
           Ptr{Complex128}, # poly
           Ptr{Cint}),# polish_only
          roots, order_changed, poly, Ref{Cint}(polish))
    roots
end

"""
    roots5(polynomial[, roots, polish=true]) -> roots

Find all the roots of `polynomial`, of degree 5 only.

Arguments:

* `polynomial`: vector of 6 coefficients (type `Number`) of the polynomial of
  which to find the roots, from the lowest coefficient to the highest one
* `roots` (optional argument): vector of initial guess roots (of length 5).  If
  you have a very rough idea where some of the roots can be, this vector is used
  as starting value for Laguerre's method
* `polish` (optional boolean keyword): if you know the roots pretty well, for
  example because you have changed the coefficiens of the polynomial only a bit,
  so the two closest roots are most likely still the closest ones, set this
  keyword to `true`.  Default it `false`

Function `roots` can be used to find roots of polynomials of any degree.
"""
function roots5{N1<:Number,N2<:Number}(poly::Vector{N1},
                                       roots::Vector{N2}=Array(Complex128,  5);
                                       polish::Bool=false)
    @assert length(poly) == 6 "Use `roots' function for polynomials of degree != 5"
    @assert length(roots) == 5 "`roots' vector must have 5 elements"
    return roots5!(float(complex(roots)), float(complex(poly)), polish)
end

end # module
