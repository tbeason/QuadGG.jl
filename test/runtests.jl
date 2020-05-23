using QuadGG
using Test

@testset "QuadGG" begin
    @testset "GaussLobatto" begin
        @testset "Continuous" begin
            @test adaptlob(sqrt,0,1;tol=1e-8) ≈ 2/3 atol=1e-8
            @test !isapprox(adaptlob(sqrt,0,1;tol=1e-3),2/3, atol=1e-8)
            @test adaptlob(x -> exp(im*x), 0,1) ≈ (exp(1im)-1)/im
        end

        @testset "Discontinuous" begin
            fun1(x) = ifelse(x<1/2,0,1)
            @test adaptlob(fun1,0,1;tol=1e-8) ≈ 1/2 atol=1e-8
        end
    end
    @testset "Simpson" begin
        @testset "Continuous" begin
            @test adaptsim(sqrt,0,1;tol=1e-8) ≈ 2/3 atol=1e-7
            @test !isapprox(adaptsim(sqrt,0,1;tol=1e-3),2/3, atol=1e-8)
            @test adaptsim(x -> exp(im*x), 0,1) ≈ (exp(1im)-1)/im
        end

        @testset "Discontinuous" begin
            fun1(x) = ifelse(x<1/2,0,1)
            @test adaptsim(fun1,0,1;tol=1e-8) ≈ 1/2 atol=1e-5

            fun2(x) = ifelse(x < 1,x+1, ifelse(x>3, 2,3-x))
            @test adaptsim(fun2,0,5;tol=1e-6,trace=true) ≈ 7.5 atol=1e-5
        end
    end
end