using LocalHilbertBurch
using Oscar
using Test

@testset "LocalHilbertBurch.jl" begin
    @test U_matrix([1,1,5]) == matrix(ZZ, [1 0 3; 1 0 3; 2 1 4; -1 -2 1])
    @test order_bounds([7]) == transpose(matrix(ZZ, [8; 1]'))
    @test degree_bounds([1,2,4]) == matrix(ZZ, [0 0 0; 0 0 0; 0 0 1; 0 0 1])
end
