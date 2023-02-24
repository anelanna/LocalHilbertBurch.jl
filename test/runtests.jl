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
    end
end
