### runtests.jl --- Test suite for PolynomialRoots.jl

# Copyright (C) 2016, 2017  Mosè Giordano

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

module PolynomialRootsTests

using PolynomialRoots
using Test

import PolynomialRoots: evalpoly

@testset "Helper functions" begin
    p, dp = @inferred(PolynomialRoots.eval_poly_der(complex(1, 2), [im, 2, 3im, 4], 3, complex(0)))
    @test p == -54 - 12im
    @test dp == -46 + 54im
    p, dp, dp2 = @inferred(PolynomialRoots.eval_poly_der2(complex(1, 2), [im, 2, 3im, 4], 3, complex(0)))
    @test p == -54 - 12im
    @test dp == -46 + 54im
    @test dp2 == 12 + 27im
    p, dp, ek = @inferred(PolynomialRoots.eval_poly_der_ek(complex(1, 2), [im, 2, 3im, 4], 3, complex(0)))
    @test p == -54 - 12im
    @test dp == -46 + 54im
    @test ek ≈ 214.10490215534654
    p, dp, dp2, ek = @inferred(PolynomialRoots.eval_poly_der2_ek(complex(1, 2), [im, 2, 3im, 4], 3, complex(0)))
    @test p == -54 - 12im
    @test dp == -46 + 54im
    @test dp2 == 12 + 27im
    @test ek ≈ 214.10490215534654

    @test PolynomialRoots.divide_poly_1(complex(5.,6.),
                                        [complex(14., -8.),
                                         complex(10., -36),
                                         complex(3., 4.)],
                                        2) == ([complex(1, 2),
                                                complex(3, 4)],
                                               complex(7, 8))

    @test [@inferred(PolynomialRoots.solve_quadratic_eq([-15.0, 8.0im, 1.0]))...] ≈ [-5im, -3im]

    @test [@inferred(PolynomialRoots.solve_cubic_eq([-6im, -(3 + 4im), 2im-2, 1.0]))...] ≈ [3, -2im, -1]
    @test [@inferred(PolynomialRoots.solve_cubic_eq(complex.([1.0, -3.0, 3.0, -1.0])))...] ≈ [1, 1, 1]

    @test isapprox(@inferred(PolynomialRoots.newton_spec([-1., 2im, 1.], 2, complex(1.), eps(1.0)))[1], -im, atol = 1e-7)
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

@testset "0th-order polynomials" begin
    poly = [1]
    res  = @inferred(roots(poly))
    @test ComplexF64[] == res

    res = @inferred(roots([0]))
    @test isnan(only(res))
end

@testset "1st-order polynomials" begin
    poly = [im, 1]
    res  = @inferred(roots(poly))
    @test zeros(length(res)) ≈ evalpoly(res, poly)
end

@testset "2nd-order polynomials" begin
    poly1 = [-15im, (5 - 3im), 1]
    res1  = @inferred(roots(poly1))
    poly2 = [0, complex(5, -0.5), 1]
    res2  = @inferred(roots(poly2))
    poly3 = [0, 0, 7]
    res3  = @inferred(roots(poly3))
    @test zeros(length(res1)) ≈ evalpoly(res1, poly1)
    @test zeros(length(res2)) ≈ evalpoly(res2, poly2)
    @test zeros(length(res1)) ≈ evalpoly(res3, poly3)
    # See examples at page 5 of http://www.cs.berkeley.edu/~wkahan/Qdrtcs.pdf
    @testset "Arbitrary precision" begin
        tol = 1e-68
        poly1 = [big"94906268.375", big"-189812534", big"94906265.625"]
        res1  = @inferred(roots(poly1))
        poly2 = [big"94906268.375", big"-189812534.75", big"94906266.375"]
        res2  = @inferred(roots(poly2))
        @test isapprox(zeros(length(res1)), evalpoly(res1, poly1), atol = tol)
        @test isapprox(zeros(length(res2)), evalpoly(res2, poly2), atol = tol)
    end
end

@testset "3rd-order polynomials" begin
    poly = [24, -(6 + 28im), (7im - 4), 1]
    @test zeros(length(poly) - 1) ≈ evalpoly(@inferred(roots(poly)), poly)
    poly = [(0.69 + 4.19im), (0.88 + 9.31im), (6.33 - 3.41im), -(9.51 - 9.91im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 3e-15)
    poly = [-(6.68 + 10.0im), (8.14 - 4.79im), -(4.51 + 0.17im), -(1.52 + 3.56im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 5e-15)
    poly = [(3.34 + 0.01im), (7.33 + 2.1im), -(2.56 + 1.09im), (6.96 + 6.18im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 6e-15)
    poly = [-(0.37 - 4.7im), -(1.77 - 5.43im), (0.94 - 6.66im), -(8.8 - 6.21im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 6e-15)
end

@testset "4th-order polynomials" begin
    tol = 2e-14
    poly = [294, -(84 + 49im), (55 + 14im), -(14 + im), 1]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = tol)
    poly = [BigInt(6), -28, 496im, 8128, -33550336im]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = tol)
    poly = [-(0.82 + 3.77im), -(0.4 - 2.11im), -(2.75 + 0.7im), (0.11 - 6.02im), (4.99 + 2.91im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = tol)
    poly = [(7.72 - 8.44im), (0.76 + 0.94im), (0.1 - 1.15im), (6.79 - 0.61im), (3.02 + 5.22im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = tol)
    poly = [(3.16 - 7.9im), (5.18 + 5.38im), (2.44 - 4.25im), (4.48 + 6.45im), -(6.27 - 3.22im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = tol)
end

@testset "5th-order polynomials" begin
    poly = [1, -5, 10, -10, 5, -1]
    @test zeros(length(poly) - 1)  ≈
        evalpoly(@inferred(roots(poly,  ones(5))),  poly)
    @test zeros(length(poly) - 1) ≈
        evalpoly(@inferred(roots5(poly, ones(5))), poly)
    poly = [-(0.76 + 0.9im), (4.96 + 0.75im), (2.45 + 0.54im), (7.35 - 2.55im), (9.18 + 3.68im), (9.17 + 8.68im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 2e-14)
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots5(poly, ones(5))), poly), atol = 2e-14)
    poly = [(2.16 - 8.7im), -(0.76 - 6.27im), (5.25 + 3.87im), (8.95 - 0.94im), (2.33 + 6.43im), (2.98 + 4.26im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 3e-14)
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots5(poly, ones(5))), poly), atol = 1e-13)
    poly = [(8.2 + 1.79im), -(8.38 + 1.57im), -(7.33 - 3.94im), (4.91 - 4.37im), (5.52 + 5.68im), -(3.9 - 7.12im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 3e-14)
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots5(poly, ones(5))), poly), atol = 5e-14)
    poly = [4.3361248929369935 + 6.315313117180402im, -8.802359600921488 - 4.624815285537038im,
            -4.258366225180867 - 6.45885917408636im, -6.136886590091395 + 6.291348712216049im,
            4.432237173050201 - 6.508641321695339im, -2.7968700254289613 + 9.101325633441487im]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 3e-14)
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots5(poly, ones(5))), poly), atol = 2e-14)
    poly = [8.045115632033657 - 5.46834836620143im, -1.571921802123125 - 3.223859508085676im,
            -0.04183883206609451 + 9.329036346463724im, -4.31270714691562 + 8.003839950544727im,
            4.802489750282724 - 1.779334725486276im, -1.1896890817078862 + 9.342349550061286im]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 6e-14)
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots5(poly, ones(5))), poly), atol = 1e-13)
    tol  = 2e-12
    poly = [120im, -(184 + 90im), (138 - 57im), (54im - 15), -(6 + 9im), 1]
    res  = @inferred(roots(poly,  [im, 2, 3im, 4, 5im]))
    res5 = @inferred(roots5(poly, [im, 2, 3im, 4, 5im]))
    @test isapprox(zeros(length(res)),  evalpoly(res,  poly), atol = tol)
    @test isapprox(zeros(length(res5)), evalpoly(res5, poly), atol = tol)
end

@testset "6th-order polynomials" begin
    poly = [(2.8 + 3.8im), -(2.0 - 0.5im), -(0.8 - 4.7im), -(8.3 - 4.5im), -(5.7 - 8.9im), -(9.1 + 2.1im), -(8.2 - 8.9im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 2e-14)
    poly = [-(1.2 - 3.5im), (2.3 + 0.3im), (5.5 - 3.9im), (0.4 + 0.3im), (6.4 - 8.9im), -(0.4 - 8.4im), (7.8 - 0.2im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 4e-14)
    poly = [-(7.1 - 6.0im), (9.2 - 0.7im), -(3.3 - 7.9im), (0.3 - 0.9im), (3.3 - 2.1im), (5.0 + 3.7im), -(3.5 + 6.7im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 4e-14)
    poly = [-(6.3 + 2.6im), -(0.3 + 7.7im), -(4.9 + 5.2im), (3.0 - 4.5im), (2.5 + 2.2im), -(2.9 - 2.1im), -(3.7 - 0.7im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 4e-14)
    poly = [-(8.4 + 2.0im), (7.0 - 2.5im), (9.9 + 5.0im), -(4.1 - 9.5im), 6.0im, -(4.9 - 6.9im), (9.4 + 4.9im)]
    @test isapprox(zeros(length(poly) - 1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 4e-14)
end

@testset "9th-order polynomials" begin
    poly = big.([-648, 3132, -6534, 7737, -5744, 2779, -878, 175, -20, 1])
    @test isapprox(zeros(length(poly)-1), evalpoly(roots(poly), poly), atol = 1e-49)
    poly = [-1, -(9 + 6im), (7 + 3im), -(7 + 5im), -(8 - 10im), (1 - 5im), -(2 - 5im), -(7 - 7im), +(7 + 6im), (7 - 10im)]
    @test isapprox(zeros(length(poly)-1), evalpoly(roots(poly), poly), atol = 2e-13)
    poly = [-(3 - 2im), (4 + 4im), (4 + 5im), (6 + 10im), (3 + 4im), -(5 + 9im), -(10 - 3im), -(3 - 8im), (6 + 3im), (6 - 2im)]
    @test isapprox(zeros(length(poly)-1), evalpoly(roots(poly), poly), atol = 2e-13)
    poly = [(2 - 10im), (8 - 1im), (3 + 8im), -(10 + 9im), (5 + 5im), (9 + 4im), (9 - 10im), (4 + 4im), 8, (10 + 9im)]
    @test isapprox(zeros(length(poly)-1), evalpoly(roots(poly), poly), atol = 2e-13)
    poly = [(4 - 5im), -(7 + 6im), -(4 - 8im), -(3 - 5im), (9 - 2im), (9 - 6im), (6 - 4im), -(8 - 6im), -(3 - 6im), (1 - 10im)]
    @test isapprox(zeros(length(poly)-1), evalpoly(roots(poly), poly), atol = 2e-13)
end

@testset "Bug fixes" begin
    # Polynomial reported in https://github.com/giordano/PolynomialRoots.jl/issues/2
    poly = [63.22256105356723, 271.8182248226112, 144.3588342039991, -25.60790629850817,
            952.6388846106129, -32.65159275777219, 62.19331611327388, 44.70637211946786,
            -28.265078307398895, -37.63653902029289, -72.26102355751738, -25.501990478720046,
            -47.40236121905153, -71.92379520244637, -57.977452001749555]
    @test isapprox(zeros(length(poly)-1),
                   evalpoly(@inferred(roots(poly)), poly), atol = 2e-11)
    # https://github.com/giordano/PolynomialRoots.jl/issues/11
    poly = [1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test isapprox(@inferred(roots(poly)), [1, 1])
    # https://github.com/giordano/PolynomialRoots.jl/pull/20
    for T in [Float32, Float64]
        poly = T[0.7513126327861701, 0.6448833539420931, 0.07782644396003469, 0.8481854810000327]
        @test isapprox(zeros(length(poly)-1), evalpoly(roots(poly), poly), atol = 1000eps(T))
    end
end

@testset "Errors" begin
    @test_throws AssertionError @inferred(roots([1,2,3],[1]))
    @test_throws AssertionError @inferred(roots5([1,2,3]))
    @test_throws AssertionError @inferred(roots5([1,2,3,4,5,6], [1]))
end

end # module
