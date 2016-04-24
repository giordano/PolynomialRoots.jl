### runtests.jl --- Test suite for CmplxRoots.jl

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

using CmplxRoots
using Base.Test

# 0th-order polynomial
poly = [1]
res  = roots(poly)
@test Array(Complex128, 0) == res

# 1st-order polynomial
tol  = 1e-12
poly = [im, 1]
res  = roots(poly)
for i = 1:length(res)
    @test_approx_eq_eps 0 (@evalpoly res[i] poly[1] poly[2]) tol
end

# 2nd-order polynomial
tol = 1e-13
poly = [-15im, (5 - 3im), 1]
res  = roots(poly)
for i = 1:length(res)
    @test_approx_eq_eps 0 (@evalpoly res[i] poly[1] poly[2] poly[3]) tol
end

# 3rd-order polynomial
tol  = 1e-13
poly = [24, -(6 + 28im), (7im - 4), 1]
res  = roots(poly)
for i = 1:length(res)
    @test_approx_eq_eps 0 (@evalpoly res[i] poly[1] poly[2] poly[3] poly[4]) tol
end

# 4th-order polynomial
tol = 1e-13
poly = [294, -(84 + 49im), (55 + 14im), -(14 + im), 1]
res  = roots(poly)
for i = 1:length(res)
    @test_approx_eq_eps 0 (@evalpoly res[i] poly[1] poly[2] poly[3] poly[4] poly[5]) tol
end

# First 5th-order polynomial
poly = [1, -5, 10, -10, 5, -1]
res  = roots(poly,  ones(5))
res5 = roots5(poly, ones(5))
for i = 1:length(res)
    @test_approx_eq 0 (@evalpoly res[i]  poly[1] poly[2] poly[3] poly[4] poly[5] poly[6])
    @test_approx_eq 0 (@evalpoly res5[i] poly[1] poly[2] poly[3] poly[4] poly[5] poly[6])
end

# Second 5th-order polynomial
tol  = 2e-12
poly = [120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1]
res  = roots(poly,  [im, 2, 3im, 4, 5im])
res5 = roots5(poly, [im, 2, 3im, 4, 5im])
for i = 1:length(res)
    @test_approx_eq_eps 0 (@evalpoly res[i]  poly[1] poly[2] poly[3] poly[4] poly[5] poly[6]) tol
    @test_approx_eq_eps 0 (@evalpoly res5[i] poly[1] poly[2] poly[3] poly[4] poly[5] poly[6]) tol
end

# Random 5th-order polynomial
info("Testing random polynomial...")
tol  = 1e-11
poly = rand(Complex128, 6)*20 - complex(10, 10)
for i = 1:length(poly); println(" poly[$i] = ", poly[i]); end
res  = roots(poly)
res5 = roots5(poly)
for i = 1:length(res)
    @test_approx_eq_eps 0 (@evalpoly res[1]  poly[1] poly[2] poly[3] poly[4] poly[5] poly[6]) tol
    @test_approx_eq_eps 0 (@evalpoly res5[1] poly[1] poly[2] poly[3] poly[4] poly[5] poly[6]) tol
end

# Throw errors
@test_throws AssertionError roots([1,2,3],[1])
@test_throws AssertionError roots5([1,2,3])
@test_throws AssertionError roots5([1,2,3,4,5,6], [1])

### Test helper functions
@test CmplxRoots.divide_poly_1(complex(5.,6.),
                               [complex(14., -8.),
                                complex(10., -36),
                                complex(3., 4.)],
                               2) == ([complex(1, 2),
                                       complex(3, 4)],
                                      complex(7, 8))

@test CmplxRoots.solve_quadratic_eq([-15., 8.im, 1.]) == (-5im, -3im)

x1, x2, x3 = CmplxRoots.solve_cubic_eq([-6im, -(3 + 4im), 2im-2, 1.])
@test_approx_eq x1  3
@test_approx_eq x2 -2im
@test_approx_eq x3 -1

@test_approx_eq_eps CmplxRoots.cmplx_newton_spec([-1., 2im, 1.], 2, complex(1.))[1] -im 1e-7
@test CmplxRoots.cmplx_newton_spec(complex([6., -5., 1.]), 2, complex(2.8))[1] == 3

@test CmplxRoots.cmplx_laguerre([-1., 2im, 1.], 2, complex(1.))[1] == -im
@test CmplxRoots.cmplx_laguerre(complex([6., -5., 1.]), 2, complex(2.8))[1] == 3

@test CmplxRoots.find_2_closest_from_5(complex([1.,3,6,10,15])) == (1,2,4.0)
@test CmplxRoots.find_2_closest_from_5(complex([1.,3,5,7,9]))   == (4,5,4.0)

a = complex([18., 5., 7., 10., 1.])
@test CmplxRoots.sort_5_points_by_separation_i(a) == [1, 5, 4, 2, 3]

@test CmplxRoots.sort_5_points_by_separation!(a) ==
    complex([18., 1., 10., 5., 7.])

@test CmplxRoots.cmplx_laguerre2newton([-1., 2im, 1.], 2, complex(1.), 2)[1] == -im
@test CmplxRoots.cmplx_laguerre2newton(complex([6., -5., 1.]), 2, complex(2.8), 2)[1] == 3
