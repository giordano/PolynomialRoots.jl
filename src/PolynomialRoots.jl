### PolynomialRoots.jl --- Fast Complex Polynomial Root Finder
#
# Copyright (C) 2012: Jan Skowron & Andrew Gould
# Copyright (C) 2016: Mosè Giordano
#
# Author: Mosè Giordano <mose AT gnu DOT org>
# Keywords: polynomials, root finding, multiple precision
#
### Licenses
#
# PolynomialRoots.jl is released under terms of Apache License, Version 2.0 or
# GNU Lesser General Public License, at your option.  Here follow license
# notices.
#
## Apache 2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
## LGPLv3+
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
### Commentary:
#
# This program started as a translation into Julia of the Fortran package
# "General Complex Polynomial Root Solver"
# (http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/) by Jan Skowron and
# Andrew Gould.
#
# A custom in the scientific comunity is (regardless of the licence you chose to
# use or distribute this software under) that if this code was important in the
# scientific process or for the results of your scientific work, we kindly ask
# you for the appropriate citation of the paper Skowron & Gould 2012 (see
# below), and we would be greatful if you pass the information about the proper
# citation to anyone whom you redistribute this software to.
#
### References:
#
#  * Skowron, J. & Gould, A. 2012, "General Complex Polynomial Root Solver and
#    Its Further Optimization for Binary Microlenses", arXiv:1203.1034.
#    URL: http://arxiv.org/abs/1203.1034
#  * Adams, D. A., 1967, "A stopping criterion for polynomial root finding",
#    Communications of the ACM, Volume 10, Issue 10, Oct. 1967, p. 655
#    DOI:10.1145/363717.363775.
#    URL: http://i.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
#
### Code:

__precompile__()

module PolynomialRoots

export roots, roots5

const third = 1//3
const MAX_ITERS = 50
const FRAC_JUMP_EVERY = 10
const FRAC_JUMPS = [0.64109297, # some random numbers
                    0.91577881, 0.25921289,  0.50487203,
                    0.08177045, 0.13653241,  0.306162  ,
                    0.37794326, 0.04618805,  0.75132137]
const FRAC_JUMP_LEN = length(FRAC_JUMPS)

function divide_poly_1{T<:AbstractFloat}(p::Complex{T},
                                         poly::Vector{Complex{T}},
                                         degree::Integer)
    coef = poly[degree+1]
    polyout = poly[1:degree]
    @inbounds for i = degree:-1:1
        prev = polyout[i]
        polyout[i] = coef
        coef = prev + p*coef
    end
    remainder = coef
    return polyout, remainder
end

function solve_quadratic_eq{T<:AbstractFloat}(poly::Vector{Complex{T}})
    c, b, a = poly
    Δ = sqrt(b*b - 4*a*c)
    if real(conj(b)*Δ) >= 0
        x0 = -0.5*(b + Δ)
    else
        x0 = -0.5*(b - Δ)
    end
    if x0 == 0
        x1 = x0
    else
        x1 = c / x0 # Viete's formula
        x0 = x0 / a
    end
    return x0, x1
end

function solve_cubic_eq{T<:AbstractFloat}(poly::Vector{Complex{T}})
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    a1  =  1 / poly[4]
    E1  = -poly[3]*a1
    E2  =  poly[2]*a1
    E3  = -poly[1]*a1
    s0  =  E1
    E12 =  E1*E1
    A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
    B   =  E12 - 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = sqrt(A*A - 4*B*B*B)
    if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s1 = (0.5*(A + Δ))^(third)
    else
        s1 = (0.5*(A - Δ))^(third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2)
end

function newton_spec{T<:AbstractFloat,E<:AbstractFloat}(poly::Vector{Complex{T}},
                                                        degree::Integer,
                                                        root::Complex{T},
                                                        epsilon::E)
    iter = 0
    success = true
    good_to_go = false
    stopping_crit2 = zero(promote_type(T,E))
    for i = 1:MAX_ITERS
        # Prepare stoping criterion.  Calculate value of polynomial and its
        # first two derivatives
        p  = poly[degree+1]
        dp = zero(Complex{T})
        if mod(i, 10) == 1 # Calculate stopping criterion every ten iterations
            ek = abs(poly[degree + 1])
            absroot = abs(root)
            # Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to
            # evaluate polynomials and derivatives
            @inbounds for k = degree:-1:1
                dp = p + dp*root
                p  = poly[k] + p*root # b_k
                # Adams (1967), equation (8).
                ek = absroot*ek + abs(p)
            end
            stopping_crit2 = abs2(epsilon*ek)
        else # Calculate just the value and derivative
            # Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to
            # evaluate polynomials and derivatives
            @inbounds for k = degree:-1:1
                dp = p + dp*root
                p  = poly[k] + p*root # b_k
            end
        end
        iter += 1
        abs2p = abs2(p)
        if abs2p == 0
            return root, iter, success
        elseif abs2p < stopping_crit2 # Simplified a little eq (10) of Adams (1967)
            if dp == 0
                return root, iter, success # if we have problem with zero, but
                                           # we are close to the root, just
                                           # accept
            end
            # do additional iteration if we are less than 10x from stopping criterion
            if abs2p < stopping_crit2*0.01
                return root, iter, success # return immediately, because we are at very good place
            else
                good_to_go = true # do one iteration more
            end
        else
            good_to_go = false # reset if we are outside the zone of the root
        end # if abs2p == 0
        if dp == 0
            # problem with zero.  Make some random jump
            dx::Complex{T} = (abs(root) + 1) * cis(FRAC_JUMPS[trunc(Integer, mod(i, FRAC_JUMP_LEN)) + 1]*2*pi)
        else
            dx = p / dp # Newton method, see http://en.wikipedia.org/wiki/Newton's_method
        end
        newroot = root - dx
        if newroot == root
            return root, iter, success # nothing changes -> return
        end
        if good_to_go # this was jump already after stopping criterion was met
            root = newroot
            return root, iter, success
        end
        if mod(i, FRAC_JUMP_EVERY) == 0 # Decide whether to do a jump of
                                        # modified length (to break cycles)
            faq = FRAC_JUMPS[trunc(Integer, mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN)) + 1]
            newroot = root - faq*dx # do jump of some semi-random length (0<faq<1)
        end
        root = newroot
    end # for i
    success = false
    return root, iter, success
end

function laguerre{T<:AbstractFloat,E<:AbstractFloat}(poly::Vector{Complex{T}},
                                                     degree::Integer,
                                                     root::Complex{T}, epsilon::E)
    c_zero = zero(Complex{T})
    iter = 0
    success = true
    good_to_go = false
    one_nth = 1 / degree
    n_1_nth = (degree - 1)*one_nth
    two_n_div_n_1 = 2 / n_1_nth
    c_one_nth = complex(one_nth)
    for i = 1:MAX_ITERS
        # prepare stoping criterion
        ek = abs(poly[degree + 1])
        absroot = abs(root)
        # calculate value of polynomial and its first two derivatives
        p = poly[degree + 1]
        dp = c_zero
        d2p_half = c_zero
        @inbounds for k = degree:-1:1 # Horner Scheme, see for eg.  Numerical
                                      # Recipes Sec. 5.3 how to evaluate
                                      # polynomials and derivatives
            d2p_half = dp + d2p_half*root
            dp = p + dp*root
            p  = poly[k] + p*root # b_k
            # Adams (1967), equation (8).
            ek = absroot*ek + abs(p)
        end
        iter=iter+1
        abs2p = abs2(p)
        if abs2p == 0
            return root, iter, success
        end
        stopping_crit2 = abs2(epsilon*ek)
        if abs2p < stopping_crit2 # simplified a little Eq. 10 of Adams (1967)
            # do additional iteration if we are less than 10x from stopping criterion
            if abs2p < 0.01*stopping_crit2
                return root, iter, success # return immediately, because we are at very good place
            else
                good_to_go= true # do one iteration more
            end
        else
            good_to_go = false # reset if we are outside the zone of the root
        end
        denom = c_zero
        if dp != 0
            fac_netwon = p / dp
            fac_extra = d2p_half / dp
            F_half = fac_netwon*fac_extra
            denom_sqrt = sqrt(1 - two_n_div_n_1*F_half)
            if real(denom_sqrt) >= 0
                denom = c_one_nth + n_1_nth*denom_sqrt
            else
                denom = c_one_nth - n_1_nth*denom_sqrt
            end
        end
        if denom == 0  # test if demoninators are > 0.0 not to divide by zero
            dx::Complex{T} = (absroot + 1) * cis(FRAC_JUMPS[trunc(Integer, mod(i,FRAC_JUMP_LEN)) + 1] * 2 * pi) # make some random jump
        else
            dx = fac_netwon / denom
        end
        newroot = root - dx
        if newroot==root
            return root, iter, success # nothing changes -> return
        end
        if good_to_go # this was jump already after stopping criterion was met
            root = newroot
            return root, iter, success
        end
        if mod(i, FRAC_JUMP_EVERY) == 0 # decide whether to do a jump of modified length (to break cycles)
            faq = FRAC_JUMPS[trunc(Integer, mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN)) + 1]
            newroot = root - faq*dx # do jump of some semi-random length (0<faq<1)
        end
        root = newroot
    end # for k
    success = false
    return root, iter, success
end

function laguerre2newton{T<:AbstractFloat,E<:AbstractFloat}(poly::Vector{Complex{T}},
                                                            degree::Integer,
                                                            root::Complex{T},
                                                            starting_mode::Integer,
                                                            epsilon::E)
    c_zero = zero(Complex{T})
    c_one = one(Complex{T})
    iter=0
    success = true
    stopping_crit2 = zero(promote_type(T, E))
    i = 0
    j = 1
    good_to_go = false
    mode = starting_mode  # mode=2 full laguerre, mode=1 SG, mode=0 newton
    # infinite loop, just to be able to come back from newton, if more than 10
    # iteration there
    while true
        #------------------------------------------------------------- mode 2
        if mode >= 2 # LAGUERRE'S METHOD
            one_nth = 1 / degree
            n_1_nth = (degree - 1)*one_nth
            two_n_div_n_1 = 2 / n_1_nth
            c_one_nth = complex(one_nth)
            for i = 1:MAX_ITERS
                # prepare stoping criterion
                ek = abs(poly[degree + 1])
                absroot = abs(root)
                # calculate value of polynomial and its first two derivatives
                p = poly[degree + 1]
                dp = c_zero
                d2p_half = c_zero
                @inbounds for k = degree:-1:1 # Horner Scheme, see for eg.
                                              # Numerical Recipes Sec. 5.3 how
                                              # to evaluate polynomials and
                                              # derivatives
                    d2p_half = dp + d2p_half*root
                    dp = p + dp*root
                    p  = poly[k] + p*root # b_k
                    # Adams (1967), equation (8).
                    ek = absroot*ek + abs(p)
                end
                abs2p = abs2(p)
                iter = iter + 1
                if abs2p == 0
                    return root, iter, success
                end
                stopping_crit2 = abs2(epsilon*ek)
                if abs2p < stopping_crit2 # (simplified a little Eq. 10 of Adams 1967)
                    # do additional iteration if we are less than 10x from stopping criterion
                    if abs2p < 0.01stopping_crit2 # ten times better than stopping criterion
                        return root, iter, success # return immediately, because we are at very good place
                    else
                        good_to_go = true # do one iteration more
                    end
                else
                    good_to_go = false # reset if we are outside the zone of the root
                end
                denom = c_zero
                if dp != 0
                    fac_netwon = p / dp
                    fac_extra = d2p_half / dp
                    F_half = fac_netwon*fac_extra
                    abs2_F_half = abs2(F_half)
                    if abs2_F_half <= 0.0625 # F<0.50, F/2<0.25
                        # go to SG method
                        if abs2_F_half <= 0.000625 # F<0.05, F/2<0.025
                            mode = 0 # go to Newton's
                        else
                            mode = 1 # go to SG
                        end
                    end
                    denom_sqrt = sqrt(c_one - two_n_div_n_1*F_half)
                    if real(denom_sqrt) >= 0
                        denom = c_one_nth + n_1_nth*denom_sqrt
                    else
                        denom = c_one_nth - n_1_nth*denom_sqrt
                    end
                end
                if denom == 0 #test if demoninators are > 0.0 not to divide by zero
                    dx = (abs(root) + 1) * cis(FRAC_JUMPS[trunc(Integer, mod(i,FRAC_JUMP_LEN)) + 1] * 2 * pi) # make some random jump
                else
                    dx = fac_netwon / denom
                end
                newroot = root - dx
                if newroot == root
                    return root, iter, success # nothing changes -> return
                end
                if good_to_go # this was jump already after stopping criterion was met
                    root = newroot
                    return root, iter, success
                end
                if mode !=2
                    root = newroot
                    j=i+1 # remember iteration index
                    break # go to Newton's or SG
                end
                if mod(i, FRAC_JUMP_EVERY) == 0 # decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS[trunc(Integer, mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN)) + 1]
                    newroot = root - faq*dx # do jump of some semi-random length (0<faq<1)
                end
                root = newroot
            end # do mode 2
            if i >= MAX_ITERS
                success = false
                return root, iter, success
            end
        end # if mode 2
        #------------------------------------------------------------- mode 1
        if mode == 1 # SECOND-ORDER GENERAL METHOD (SG)
            for i = j:MAX_ITERS
                # calculate value of polynomial and its first two derivatives
                p = poly[degree + 1]
                dp = c_zero
                d2p_half = c_zero
                if mod(i - j, 10) == 0
                    # prepare stoping criterion
                    ek = abs(poly[degree+1])
                    absroot = abs(root)
                    @inbounds for k = degree:-1:1 # Horner Scheme, see for eg.
                                                  # Numerical Recipes Sec. 5.3
                                                  # how to evaluate polynomials
                                                  # and derivatives
                        d2p_half = dp + d2p_half*root
                        dp = p + dp*root
                        p  = poly[k] + p*root # b_k
                        # Adams (1967) equation (8).
                        ek = absroot*ek + abs(p)
                    end
                    stopping_crit2 = abs2(epsilon*ek)
                else
                    @inbounds for k = degree:-1:1 # Horner Scheme, see for eg.
                                                  # Numerical Recipes Sec. 5.3
                                                  # how to evaluate polynomials
                                                  # and derivatives
                        d2p_half = dp + d2p_half*root
                        dp = p + dp*root
                        p  = poly[k] + p*root # b_k
                    end
                end
                abs2p = abs2(p) #abs(p)**2
                iter = iter + 1
                if abs2p == 0
                    return root, iter, success
                end
                if abs2p < stopping_crit2 # (simplified a little Eq. 10 of Adams 1967)
                    if dp == 0
                        return root, iter, success
                    end
                    # do additional iteration if we are less than 10x from stopping criterion
                    if abs2p < 0.01*stopping_crit2 # ten times better than stopping criterion
                        return root, iter, success # return immediately, because we are at very good place
                    else
                        good_to_go = true # do one iteration more
                    end
                else
                    good_to_go = false # reset if we are outside the zone of the root
                end
                if dp == 0 #test if demoninators are > 0.0 not to divide by zero
                    dx = (abs(root) + 1) * cis(FRAC_JUMPS[trunc(Integer, mod(i,FRAC_JUMP_LEN)) + 1]*2*pi) # make some random jump
                else
                    fac_netwon = p / dp
                    fac_extra = d2p_half / dp
                    F_half = fac_netwon*fac_extra
                    abs2_F_half = abs2(F_half)
                    if abs2_F_half <= 0.000625 # F<0.05, F/2<0.025
                        mode = 0 # set Newton's, go there after jump
                    end
                    dx = fac_netwon*(c_one + F_half)  # SG
                end
                newroot = root - dx
                if newroot == root
                    return root, iter, success # nothing changes -> return
                end
                if good_to_go # this was jump already after stopping criterion was met
                    root = newroot
                    return root, iter, success
                end
                if mode !=1
                    root = newroot
                    j = i + 1 # remember iteration number
                    break # go to Newton's
                end
                if mod(i, FRAC_JUMP_EVERY) == 0 # decide whether to do a jump of modified length (to break cycles)
                    faq = FRAC_JUMPS[trunc(Integer, mod(i / FRAC_JUMP_EVERY - 1, FRAC_JUMP_LEN)) + 1]
                    newroot = root - faq*dx # do jump of some semi-random length (0<faq<1)
                end
                root = newroot
            end # do mode 1
            if i >= MAX_ITERS
                success = false
                return root, iter, success
            end
        end # if mode 1
        if mode == 0 # NEWTON'S METHOD
            for i = j:j+10 # do only 10 iterations the most, then go back to full Laguerre's
                # calculate value of polynomial and its first two derivatives
                p = poly[degree + 1]
                dp = c_zero
                if i == j # calculate stopping crit only once at the begining
                    # prepare stoping criterion
                    ek = abs(poly[degree+1])
                    absroot = abs(root)
                    @inbounds for k = degree:-1:1 # Horner Scheme, see for eg.
                                                  # Numerical Recipes Sec. 5.3
                                                  # how to evaluate polynomials
                                                  # and derivatives
                        dp = p + dp*root
                        p  = poly[k] + p*root # b_k
                        # Adams (1967), equation (8).
                        ek = absroot*ek + abs(p)
                    end
                    stopping_crit2 = abs2(epsilon*ek)
                else #
                    @inbounds for k = degree:-1:1 # Horner Scheme, see for eg.
                                                  # Numerical Recipes Sec. 5.3
                                                  # how to evaluate polynomials
                                                  # and derivatives
                        dp = p + dp*root
                        p = poly[k] + p*root # b_k
                    end
                end
                abs2p = abs2(p) #abs(p)**2
                iter = iter + 1
                if abs2p == 0
                    return root, iter, success
                end
                if abs2p < stopping_crit2 # (simplified a little Eq. 10 of Adams 1967)
                    if dp == 0
                        return root, iter, success
                    end
                    # do additional iteration if we are less than 10x from stopping criterion
                    if abs2p < 0.01*stopping_crit2 # ten times better than stopping criterion
                        return root, iter, success # return immediately, because we are at very good place
                    else
                        good_to_go = true # do one iteration more
                    end
                else
                    good_to_go = false # reset if we are outside the zone of the root
                end
                if dp == 0 # test if demoninators are > 0.0 not to divide by zero
                    dx = (abs(root) + 1) * cis(FRAC_JUMPS[trunc(Integer, mod(i,FRAC_JUMP_LEN)) + 1] * 2 * pi) # make some random jump
                else
                    dx = p / dp
                end
                newroot = root - dx
                if newroot == root
                    return root, iter, success # nothing changes -> return
                end
                if good_to_go
                    root = newroot
                    return root, iter, success
                end
                root = newroot
            end # do mode 0 10 times
            if iter >= MAX_ITERS
                # too many iterations here
                success=false
                return root, iter, success
            end
            mode = 2 # go back to Laguerre's. This happens when we were unable
                     # to converge in 10 iterations with Newton's
        end # if mode 0
    end # end of infinite loop (while true)
    success = false
    return root, iter, success
end

function find_2_closest_from_5{T<:AbstractFloat}(points::Vector{Complex{T}})
    n = 5
    d2min = Inf
    i1 = 0
    i2 = 0
    for j = 1:n, i = 1:j-1
        d2 = abs2(points[i] - points[j])
        if d2 <= d2min
            i1 = i
            i2 = j
            d2min = d2
        end
    end
    return i1, i2, d2min
end

function sort_5_points_by_separation_i{T<:AbstractFloat}(points::Vector{Complex{T}})
    n = 5
    distances2 = ones(T, n, n)*Inf
    dmin = Array{T}(n)
    for j = 1:n, i = 1:j-1
        distances2[i, j] = distances2[j, i] = abs2(points[i] - points[j])
    end
    @inbounds for j = 1:n
        dmin[j] = minimum(distances2[j,:])
    end
    return sort(collect(1:n), lt=(i,j) -> dmin[i]>dmin[j])
end

function sort_5_points_by_separation!{T<:AbstractFloat}(points::Vector{Complex{T}})
    points = points[sort_5_points_by_separation_i(points)]
end

# Original function has a `use_roots_as_starting_points' argument.  We don't
# have this argument and always use `roots' as starting points, it's a task of
# the interface to set a proper starting value if the user doesn't provide it.
function roots!{T<:AbstractFloat,E<:AbstractFloat}(roots::Vector{Complex{T}},
                                                   poly::Vector{Complex{T}},
                                                   epsilon::E, degree::Integer,
                                                   polish::Bool)
    isnan(epsilon) && (epsilon = eps(T))
    poly2 = copy(poly)
    # skip small degree polynomials from doing Laguerre's method
    if degree <= 1
        if degree == 1
            roots[1] = -poly[1] / poly[2]
        end
        return roots
    end
    @inbounds for n = degree:-1:3
        roots[n], iter, success = laguerre2newton(poly2, n, roots[n], 2, epsilon)
        if ! success
            roots[n], iter, success = laguerre(poly2, n,
                                               zero(Complex{T}), epsilon)
        end
        # divide the polynomial by this root
        coef = poly2[n+1]
        @inbounds for i = n:-1:1
            prev = poly2[i]
            poly2[i] = coef
            coef = prev + roots[n]*coef
        end
    end
    # Differently from original function, we always calculate last 2 roots with
    # `solve_quadratic_eq'.
    roots[1], roots[2] = solve_quadratic_eq(poly2)
    if polish
        @inbounds for n = 1:degree # polish roots one-by-one with a full polynomial
            roots[n], iter, success = laguerre(poly, degree, roots[n], epsilon)
        end
    end
    return roots
end

function roots{N1<:Number,N2<:Number,E<:AbstractFloat}(poly::Vector{N1},
                                                       roots::Vector{N2};
                                                       epsilon::E=NaN,
                                                       polish::Bool=false)
    degree = length(poly) - 1
    @assert degree == length(roots) "`poly' must have one element more than `roots'"
    roots!(promote(float.(complex(roots)), float.(complex(poly)))...,
           epsilon, degree, polish)
end

function roots{N<:Number,E<:AbstractFloat}(poly::Vector{N};
                                           epsilon::E=NaN,
                                           polish::Bool=false)
    degree = length(poly) - 1
    roots!(zeros(Complex{real(float(N))}, degree), float.(complex(poly)),
           epsilon, degree, polish)
end

function roots5!{T<:AbstractFloat,E<:AbstractFloat}(roots::Vector{Complex{T}},
                                                    poly::Vector{Complex{T}},
                                                    epsilon::E, polish::Bool)
    isnan(epsilon) && (epsilon = eps(T))
    c_zero = zero(Complex{T})
    degree = 5
    roots_robust = copy(roots)
    go_to_robust = 0
    if ! polish
        # The roots are assumed to have been initialized by the user interface.
        go_to_robust = 1
    end
    first_3_roots_order_changed = false
    succ = false
    for loops = 1:3
        # ROBUST
        # (we do not know the roots)
        if go_to_robust > 0
            if go_to_robust > 2 # something is wrong
                return roots_robust # return not-polished roots, because polishing creates errors
            end
            poly2 = copy(poly) # copy coeffs
            @inbounds for m = degree:-1:4 # find the roots one-by-one (until 3 are left to be found)
                roots[m], iter, succ = laguerre2newton(poly2, m, roots[m],
                                                       2, epsilon)
                if ! succ
                    roots[m], iter, succ = laguerre(poly2, m, c_zero, epsilon)
                end
                # divide polynomial by this root
                poly2, remainder = divide_poly_1(roots[m], poly2, m)
            end
            # find last 3 roots with cubic euqation solver (Lagrange's method)
            roots[1], roots[2], roots[3] = solve_cubic_eq(poly2)
            # all roots found
            # sort roots - first will be most isolated, last two will be the closest
            sort_5_points_by_separation!(roots)
            # copy roots in case something will go wrong during polishing
            roots_robust = copy(roots)
            # set flag, that roots have been resorted
            first_3_roots_order_changed = true
        end  # go_to_robust>0
        # POLISH
        # (we know the roots approximately, and we guess that last two are closest)
        #---------------------
        poly2 = copy(poly) # copy coeffs
        @inbounds for m = 1:degree-2
            # polish roots with full polynomial
            roots[m], iter, succ = newton_spec(poly2, degree, roots[m], epsilon)
            if ! succ
                # go back to robust
                go_to_robust = go_to_robust + 1
                roots *= c_zero
                break
            end
        end # m = 1:degree-2
        if succ
            # comment out division and quadratic if you (POWN) polish with Newton only
            @inbounds for m = 1:degree-2
                poly2, remainder = divide_poly_1(roots[m], poly2, degree - m + 1)
            end
            # last two roots are found with quadratic equation solver
            # (this is faster and more robust, although little less accurate)
            roots[degree-1], roots[degree] = solve_quadratic_eq(poly2)
            # all roots found and polished
            # TEST ORDER
            # test closest roots if they are the same pair as given to polish
            root4,root5, d2min = find_2_closest_from_5(roots)
            if (root4 < degree - 1) || (root5 < degree - 1)
                # after polishing some of the 3 far roots become one of the 2
                # closest ones go back to robust
                if go_to_robust > 0
                    # if came from robust copy two most isolated roots as
                    # starting points for new robust
                    @inbounds for i = 1:degree-3
                        roots[degree-i+1] = roots_robust[i]
                    end
                else
                    # came from users initial guess copy some 2 roots (except
                    # the closest ones)
                    i2 = degree
                    @inbounds for i = 1:degree
                        if (i != root4) && (i != root5)
                            roots[i2] = roots[i]
                            i2 -= 1
                        end
                        if i2 <= 3
                            break # do not copy those that will be done by cubic in robust
                        end
                    end
                end
                go_to_robust = go_to_robust + 1
            else
                # root4 and root5 comes from the initial closest pair
                # most common case
                return roots
            end
        end
        #---------------------
    end # loops
    return roots
end

function roots5{N1<:Number,N2<:Number,E<:AbstractFloat}(poly::Vector{N1},
                                                        roots::Vector{N2};
                                                        epsilon::E=NaN)
    @assert length(poly) == 6 "Use `roots' function for polynomials of degree != 5"
    @assert length(roots) == 5 "`roots' vector must have 5 elements"
    return roots5!(promote(float.(complex(roots)), float.(complex(poly)))...,
                   epsilon, true)
end

function roots5{N<:Number,E<:AbstractFloat}(poly::Vector{N}; epsilon::E=NaN)
    @assert length(poly) == 6 "Use `roots' function for polynomials of degree != 5"
    return roots5!(zeros(Complex{real(float(N))},  5), float.(complex(poly)),
                   epsilon, false)
end

"""
    roots(polynomial[, roots]; polish=true, epsilon=1e-20) -> roots

Find all the roots of `polynomial`, of any degree.

Mandatory argument:

* `polynomial`: vector of coefficients (type `Number`) of the polynomial of
  which to find the roots, from the lowest coefficient to the highest one

Optional argument:

* `roots`: vector of initial guess roots.  If you have a very rough idea where
  some of the roots can be, this vector is used as starting value for Laguerre's
  method

Optional keywords:

* `polish` (`Bool`): if set to `true`, after all roots have been found by
  dividing original polynomial by each root found, all roots will be polished
  using full polynomial.  Default is `false`

* `epsilon` (`AbstractFloat`): this is used to determine the stopping criterion
  described in Skowron & Gould paper.  If not set, it defaults to machine
  precision of `polynomial` (and `roots`) argument(s).  This is *not* the
  precision with which the roots will be calculated

Function `root5` is specialized for polynomials of degree 5.
"""
roots

"""
    roots5(polynomial[, roots]; epsilon=1e-20) -> roots

Find all the roots of `polynomial`, of degree 5 only.

Mandatory argument:

* `polynomial`: vector of 6 coefficients (type `Number`) of the polynomial of
  which to find the roots, from the lowest coefficient to the highest one

Optional argument:

* `roots`: vector of initial guess roots (of length 5).  If you have a very
  rough idea where some of the roots can be, this vector is used as starting
  value for Laguerre's method and the provided roots will be only polished

Optional keyword:

* `epsilon` (`AbstractFloat`): this is used to determine the stopping criterion
  described in Skowron & Gould paper.  If not set, it defaults to machine
  precision of `polynomial` (and `roots`) argument(s).  This is *not* the
  precision with which the roots will be calculated

Function `roots` can be used to find roots of polynomials of any degree.
"""
roots5

# Use algorithm described in Knuth, TAOCP vol. 2, section 4.6.4, equation (3) to
# evaluate polynomial with coefficients u at point z.  This is the same
# algorithm used in @evalpoly macro.
function evalpoly{T<:AbstractFloat}(z::Complex{T}, u::Vector{Complex{T}})
    x, y = reim(z)
    r = x + x
    s = x*x + y*y
    degree = length(u) - 1
    a = similar(u, degree)
    b = similar(u, degree)
    a[1] = u[end]
    b[1] = u[end-1]
    for j = 2:degree
        a[j] = b[j-1] + r*a[j-1]
        b[j] = u[end-j] -s*a[j-1]
    end
    z*a[end] + b[end]
end
function evalpoly{Z<:AbstractFloat,U<:Number}(z::Complex{Z}, u::Vector{U})
    T = promote_type(Complex{Z}, U)
    return evalpoly(convert(T, z), convert(Vector{T}, u))
end
evalpoly{Z,U}(z::Vector{Complex{Z}}, u::Vector{U}) =
    map(x->evalpoly(x, u), z)

end # module
