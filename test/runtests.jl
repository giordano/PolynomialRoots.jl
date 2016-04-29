### runtests.jl --- Test suite for PolynomialRoots.jl

# Copyright (C) 2016  Mosè Giordano

# Maintainer: Mosè Giordano <mose AT gnu DOT org>

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

using PolynomialRoots
using Base.Test

# 0th-order polynomial
poly = [1]
res  = roots(poly)
@test Array(Complex128, 0) == res

# 1st-order polynomial
poly = [im, 1]
res  = roots(poly)
@test_approx_eq zeros(length(res)) PolynomialRoots.evalpoly(res, poly)

# 2nd-order polynomials
poly1 = [-15im, (5 - 3im), 1]
res1  = roots(poly1)
poly2 = [0, complex(5, -0.5), 1]
res2  = roots(poly2)
poly3 = [0, 0, 7]
res3  = roots(poly3)
@test_approx_eq zeros(length(res1)) PolynomialRoots.evalpoly(res1, poly1)
@test_approx_eq zeros(length(res2)) PolynomialRoots.evalpoly(res2, poly2)
@test_approx_eq zeros(length(res1)) PolynomialRoots.evalpoly(res3, poly3)

# Test multiple precision.  See examples at page 5 of
# http://www.cs.berkeley.edu/~wkahan/Qdrtcs.pdf
tol = 1e-68
poly1 = [big"94906268.375", big"-189812534", big"94906265.625"]
res1  = roots(poly1, epsilon=1e-70)
poly2 = [big"94906268.375", big"-189812534.75", big"94906266.375"]
res2  = roots(poly2, epsilon=1e-70)
@test_approx_eq_eps zeros(length(res1)) PolynomialRoots.evalpoly(res1, poly1) tol
@test_approx_eq_eps zeros(length(res2)) PolynomialRoots.evalpoly(res2, poly2) tol

# 3rd-order polynomial
poly = [24, -(6 + 28im), (7im - 4), 1]
res  = roots(poly)
@test_approx_eq zeros(length(res)) PolynomialRoots.evalpoly(res, poly)

# 4th-order polynomial
tol = 1e-13
poly1 = [294, -(84 + 49im), (55 + 14im), -(14 + im), 1]
res1  = roots(poly1)
poly2 = [BigInt(6), -28, 496im, 8128, -33550336im]
res2  = roots(poly2, epsilon=1e-14)
@test_approx_eq_eps zeros(length(res1)) PolynomialRoots.evalpoly(res1, poly1) tol
@test_approx_eq_eps zeros(length(res2)) PolynomialRoots.evalpoly(res2, poly2) tol

# First 5th-order polynomial
poly = [1, -5, 10, -10, 5, -1]
res  = roots(poly,  ones(5))
res5 = roots5(poly, ones(5))
@test_approx_eq zeros(length(res))  PolynomialRoots.evalpoly(res,  poly)
@test_approx_eq zeros(length(res5)) PolynomialRoots.evalpoly(res5, poly)

# Second 5th-order polynomial
tol  = 2e-12
poly = [120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1]
res  = roots(poly,  [im, 2, 3im, 4, 5im])
res5 = roots5(poly, [im, 2, 3im, 4, 5im])
@test_approx_eq_eps zeros(length(res))  PolynomialRoots.evalpoly(res,  poly) tol
@test_approx_eq_eps zeros(length(res5)) PolynomialRoots.evalpoly(res5, poly) tol

# Random 5th-order polynomials
info("Testing random polynomials...")
tol  = 2e-12
poly1 = rand(Complex128, 6)*20 - complex(10, 10)
for i = 1:length(poly1); println(" poly1[$i] = ", poly1[i]); end
res1  = roots(poly1, polish=true)
res51 = roots5(poly1)
poly2 = rand(Complex128, 6)*20 - complex(10, 10)
for i = 1:length(poly2); println(" poly2[$i] = ", poly2[i]); end
res2  = roots(promote(poly2, zeros(Complex{BigFloat}, 5))..., polish=true)
res52 = roots5(promote(poly2, zeros(Complex{BigFloat}, 5))...)
@test_approx_eq_eps zeros(length(res1))  PolynomialRoots.evalpoly(res1,  poly1) tol
@test_approx_eq_eps zeros(length(res51)) PolynomialRoots.evalpoly(res51, poly1) tol
@test_approx_eq_eps zeros(length(res2))  PolynomialRoots.evalpoly(res2,  poly2) tol
@test_approx_eq_eps zeros(length(res52)) PolynomialRoots.evalpoly(res52, poly2) tol

# Throw errors
@test_throws AssertionError roots([1,2,3],[1])
@test_throws AssertionError roots5([1,2,3])
@test_throws AssertionError roots5([1,2,3,4,5,6], [1])

### Test helper functions
@test PolynomialRoots.divide_poly_1(complex(5.,6.),
                                    [complex(14., -8.),
                                     complex(10., -36),
                                     complex(3., 4.)],
                                    2) == ([complex(1, 2),
                                            complex(3, 4)],
                                           complex(7, 8))

@test PolynomialRoots.solve_quadratic_eq([-15., 8.im, 1.]) == (-5im, -3im)

x1, x2, x3 = PolynomialRoots.solve_cubic_eq([-6im, -(3 + 4im), 2im-2, 1.])
@test_approx_eq x1  3
@test_approx_eq x2 -2im
@test_approx_eq x3 -1

@test_approx_eq_eps PolynomialRoots.newton_spec([-1., 2im, 1.], 2, complex(1.), eps(1.0))[1] -im 1e-7
@test PolynomialRoots.newton_spec(complex([6., -5., 1.]), 2, complex(2.8), eps(1.0))[1] == 3

@test PolynomialRoots.laguerre([-1., 2im, 1.], 2, complex(1.), eps(1.0))[1] == -im
@test PolynomialRoots.laguerre(complex([6., -5., 1.]), 2, complex(2.8), eps(1.0))[1] == 3

@test PolynomialRoots.find_2_closest_from_5(complex([1.,3,6,10,15])) == (1,2,4.0)
@test PolynomialRoots.find_2_closest_from_5(complex([1.,3,5,7,9]))   == (4,5,4.0)

a = complex([18., 5., 7., 10., 1.])
@test PolynomialRoots.sort_5_points_by_separation_i(a) == [1, 5, 4, 2, 3]

@test PolynomialRoots.sort_5_points_by_separation!(a) ==
    complex([18., 1., 10., 5., 7.])

@test PolynomialRoots.laguerre2newton([-1., 2im, 1.], 2, complex(1.), 2, eps(1.0))[1] == -im
@test PolynomialRoots.laguerre2newton(complex([6., -5., 1.]), 2, complex(2.8), 2, eps(1.0))[1] == 3
