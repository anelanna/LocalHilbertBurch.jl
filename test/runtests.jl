using LocalHilbertBurch
using Oscar
using Test

@testset "LocalHilbertBurch.jl" begin
    @testset "helper methods" begin
        @test U_matrix([1,1,5]) == matrix(ZZ, [1 0 3; 1 0 3; 2 1 4; -1 -2 1])
        @test order_bounds([7]) == transpose(matrix(ZZ, [8; 1]'))
        @test degree_bounds([1,2,4]) == matrix(ZZ, [0 0 0; 0 0 0; 0 0 1; 0 0 1])
    end

    @testset "cells" begin
        list = Dict{Int64, Dict{Int64, Vector{Cell}}}()
        for n in 1:5
            list[n] = sorted_celllist(n)
        end
        @test length(list) == 5
        @test [length(list[i]) for i in 1:5] == [1,2,3,4,5]
        R1,_ = PolynomialRing(QQ, [["x","y"];["c_"*string(i) for i in 1:23]])
        Q,(x,y) = GradedPolynomialRing(QQ, ["x", "y"])
        C1 = Cell([1,5,8,10],R1,Q)
        @test C1.hilb == [1,2,3,4,3,3,3,2,2,1]
        @test C1.U == matrix(ZZ,[1    4   6   7; 1    4   6   7;-2    1   3   4;-4   -1   1   2;-5   -2   0   1])
        @test C1.H == matrix(Q,[y     0     0     0;-x   y^4     0     0; 0    -x   y^3     0; 0     0    -x   y^2; 0     0     0    -x])
        @test C1.M == matrix(R1,[y 0 0 0;  -x y^4 0 0; c_1 -x + y^3*c_4 + y^2*c_3 + y*c_2   y^3 0; c_5 y^3*c_9 + y^2*c_8 + y*c_7 + c_6     -x + y^2*c_11 + y*c_10           y^2;c_12   y^3*c_16 + y^2*c_15 + y*c_14 + c_13   y^2*c_19 + y*c_18 + c_17   -x + y*c_20])
        R2,_ = PolynomialRing(QQ, [["x","y"];["c_"*string(i) for i in 1:17]])
        C2 = Cell([2,3,5,7],R2,Q)
        @test C2.hilb == [1, 2, 3, 4, 4, 2, 1]
        @test C2.H == matrix(Q,[y^2    0     0     0; -x    y     0     0;  0   -x   y^2     0; 0    0    -x   y^2;  0    0     0    -x])
        @test C2.U == matrix(ZZ, [ 2    2   3   4; 1    1   2   3; 1    1   2   3; 0    0   1   2;-1   -1   0   1])
        @test C2.M == matrix(R2,[        y^2     0               0             0; -x + y*c_1     y               0             0;      y*c_2    -x             y^2             0;y*c_4 + c_3   c_5      -x + y*c_6           y^2;y*c_8 + c_7   c_9   y*c_11 + c_10   -x + y*c_12         ])
    end
end
