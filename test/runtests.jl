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
        R1,(x,y,c...) = PolynomialRing(QQ, [["x","y"];["c["*string(i)*"]" for i in 1:23]])
        Q,(x,y) = GradedPolynomialRing(QQ, ["x", "y"])
        C1 = Cell([1,5,8,10],R1,Q)
        @test C1.hilb == [1,2,3,4,3,3,3,2,2,1]
        @test C1.U == matrix(ZZ,[1    4   6   7; 1    4   6   7;-2    1   3   4;-4   -1   1   2;-5   -2   0   1])
        @test C1.H == matrix(Q,[y     0     0     0;-x   y^4     0     0; 0    -x   y^3     0; 0     0    -x   y^2; 0     0     0    -x])
        @test C1.M == matrix(R1, [y 0 0  0;   -x y^4 0 0; c[1] -x + y^3*c[4] + y^2*c[3] + y*c[2] y^3 0; c[5] y^3*c[9] + y^2*c[8] + y*c[7] + c[6] -x + y^2*c[11] + y*c[10] y^2; c[12]   y^3*c[16] + y^2*c[15] + y*c[14] + c[13]   y^2*c[19] + y*c[18] + c[17]   -x + y*c[20]])
        R2,(x,y,c...) = PolynomialRing(QQ, [["x","y"];["c["*string(i)*"]" for i in 1:17]])
        C2 = Cell([2,3,5,7],R2,Q)
        @test C2.hilb == [1, 2, 3, 4, 4, 2, 1]
        @test C2.H == matrix(Q,[y^2    0     0     0; -x    y     0     0;  0   -x   y^2     0;  0    0    -x   y^2; 0    0     0    -x])
        @test C2.U == matrix(ZZ, [ 2    2   3   4; 1    1   2   3; 1    1   2   3; 0    0   1   2;-1   -1   0   1])
        @test C2.M == matrix(R2, [y^2 0 0 0; -x + y*c[1] y 0 0; y*c[2] -x y^2 0;y*c[4] + c[3] c[5] -x + y*c[6]   y^2; y*c[8] + c[7]   c[9]   y*c[11] + c[10]   -x + y*c[12]])
    end
end
