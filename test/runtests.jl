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

@testset "0th-order polynomials" begin
    poly = [1]
    res  = @inferred(roots(poly))
    @test Array{Complex128}(0) == res
end

@testset "1st-order polynomials" begin
    poly = [im, 1]
    res  = @inferred(roots(poly))
    @test zeros(length(res)) ≈ PolynomialRoots.evalpoly(res, poly)
end

@testset "2nd-order polynomials" begin
    poly1 = [-15im, (5 - 3im), 1]
    res1  = @inferred(roots(poly1))
    poly2 = [0, complex(5, -0.5), 1]
    res2  = @inferred(roots(poly2))
    poly3 = [0, 0, 7]
    res3  = @inferred(roots(poly3))
    @test zeros(length(res1)) ≈ PolynomialRoots.evalpoly(res1, poly1)
    @test zeros(length(res2)) ≈ PolynomialRoots.evalpoly(res2, poly2)
    @test zeros(length(res1)) ≈ PolynomialRoots.evalpoly(res3, poly3)
    # See examples at page 5 of http://www.cs.berkeley.edu/~wkahan/Qdrtcs.pdf
    @testset "Arbitrary precision" begin
        tol = 1e-68
        poly1 = [big"94906268.375", big"-189812534", big"94906265.625"]
        res1  = @inferred(roots(poly1))
        poly2 = [big"94906268.375", big"-189812534.75", big"94906266.375"]
        res2  = @inferred(roots(poly2))
        @test zeros(length(res1)) ≈ PolynomialRoots.evalpoly(res1, poly1) atol = tol
        @test zeros(length(res2)) ≈ PolynomialRoots.evalpoly(res2, poly2) atol = tol
    end
end

@testset "3rd-order polynomials" begin
    poly = [24, -(6 + 28im), (7im - 4), 1]
    res  = @inferred(roots(poly))
    @test zeros(length(res)) ≈ PolynomialRoots.evalpoly(res, poly)
end

@testset "4th-order polynomials" begin
    tol = 1e-13
    poly1 = [294, -(84 + 49im), (55 + 14im), -(14 + im), 1]
    res1  = @inferred(roots(poly1))
    poly2 = [BigInt(6), -28, 496im, 8128, -33550336im]
    res2  = @inferred(roots(poly2, epsilon=1e-14))
    @test zeros(length(res1)) ≈ PolynomialRoots.evalpoly(res1, poly1) atol = tol
    @test zeros(length(res2)) ≈ PolynomialRoots.evalpoly(res2, poly2) atol = tol
end

@testset "5th-order polynomials" begin
    poly = [1, -5, 10, -10, 5, -1]
    res  = @inferred(roots(poly,  ones(5)))
    res5 = @inferred(roots5(poly, ones(5)))
    @test zeros(length(res))  ≈ PolynomialRoots.evalpoly(res,  poly)
    @test zeros(length(res5)) ≈ PolynomialRoots.evalpoly(res5, poly)

    tol  = 2e-12
    poly = [120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1]
    res  = @inferred(roots(poly,  [im, 2, 3im, 4, 5im]))
    res5 = @inferred(roots5(poly, [im, 2, 3im, 4, 5im]))
    @test zeros(length(res))  ≈ PolynomialRoots.evalpoly(res,  poly) atol = tol
    @test zeros(length(res5)) ≈ PolynomialRoots.evalpoly(res5, poly) atol = tol
end

@testset "9th-order polynomials" begin
    poly = big.([-648, 3132, -6534, 7737, -5744, 2779, -878, 175, -20, 1])
    @test zeros(length(poly)-1) ≈ PolynomialRoots.evalpoly(roots(poly), poly) atol = 1e-49
end

@testset "Bug fixes" begin
    # Polynomial reported in https://github.com/giordano/PolynomialRoots.jl/issues/2
    poly = [63.22256105356723, 271.8182248226112, 144.3588342039991, -25.60790629850817,
            952.6388846106129, -32.65159275777219, 62.19331611327388, 44.70637211946786,
            -28.265078307398895, -37.63653902029289, -72.26102355751738, -25.501990478720046,
            -47.40236121905153, -71.92379520244637, -57.977452001749555]
    @test zeros(length(poly)-1) ≈
        PolynomialRoots.evalpoly(@inferred(roots(poly)), poly) atol = 2e-11
end

@testset "Random polynomials" begin
    tol  = 2e-12
    poly1 = rand(Complex128, 6)*20 - complex(10, 10)
    for i = 1:length(poly1); println(" poly1[$i] = ", poly1[i]); end
    res1  = @inferred(roots(poly1, polish=true))
    res51 = roots5(poly1)
    poly2 = rand(Complex128, 6)*20 - complex(10, 10)
    for i = 1:length(poly2); println(" poly2[$i] = ", poly2[i]); end
    res2  = @inferred(roots(promote(poly2, zeros(Complex{BigFloat}, 5))..., polish=true))
    res52 = @inferred(roots5(promote(poly2, zeros(Complex{BigFloat}, 5))...))
    @test zeros(length(res1))  ≈ PolynomialRoots.evalpoly(res1,  poly1) atol = tol
    @test zeros(length(res51)) ≈ PolynomialRoots.evalpoly(res51, poly1) atol = tol
    @test zeros(length(res2))  ≈ PolynomialRoots.evalpoly(res2,  poly2) atol = tol
    @test zeros(length(res52)) ≈ PolynomialRoots.evalpoly(res52, poly2) atol = tol
end

@testset "Errors" begin
    @test_throws AssertionError @inferred(roots([1,2,3],[1]))
    @test_throws AssertionError @inferred(roots5([1,2,3]))
    @test_throws AssertionError @inferred(roots5([1,2,3,4,5,6], [1]))
end

@testset "Helper functions" begin
    @test PolynomialRoots.divide_poly_1(complex(5.,6.),
                                        [complex(14., -8.),
                                         complex(10., -36),
                                         complex(3., 4.)],
                                        2) == ([complex(1, 2),
                                                complex(3, 4)],
                                               complex(7, 8))

    @test @inferred(PolynomialRoots.solve_quadratic_eq([-15.0, 8.0im, 1.0])) == (-5im, -3im)

    x1, x2, x3 = @inferred(PolynomialRoots.solve_cubic_eq([-6im, -(3 + 4im), 2im-2, 1.]))
    @test x1 ≈  3
    @test x2 ≈ -2im
    @test x3 ≈ -1

    @test @inferred(PolynomialRoots.newton_spec([-1., 2im, 1.], 2, complex(1.), eps(1.0)))[1] ≈ -im atol = 1e-7
    @test @inferred(PolynomialRoots.newton_spec(complex([6., -5., 1.]), 2, complex(2.8), eps(1.0)))[1] == 3

    @test @inferred(PolynomialRoots.laguerre([-1., 2im, 1.], 2, complex(1.), eps(1.0)))[1] == -im
    @test @inferred(PolynomialRoots.laguerre(complex([6., -5., 1.]), 2, complex(2.8), eps(1.0)))[1] == 3

    @test @inferred(PolynomialRoots.find_2_closest_from_5(complex([1.,3,6,10,15]))) == (1,2,4.0)
    @test @inferred(PolynomialRoots.find_2_closest_from_5(complex([1.,3,5,7,9])))   == (4,5,4.0)

    a = complex([18., 5., 7., 10., 1.])
    @test @inferred(PolynomialRoots.sort_5_points_by_separation_i(a)) == [1, 5, 4, 2, 3]

    @test @inferred(PolynomialRoots.sort_5_points_by_separation!(a)) ==
        complex([18., 1., 10., 5., 7.])

    @test @inferred(PolynomialRoots.laguerre2newton([-1., 2im, 1.], 2, complex(1.), 2, eps(1.0)))[1] == -im
    @test @inferred(PolynomialRoots.laguerre2newton(complex([6., -5., 1.]), 2, complex(2.8), 2, eps(1.0)))[1] == 3
end
