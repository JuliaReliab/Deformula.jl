import Deformula: _deint, deformulaZeroToInf, deformulaMinusOneToOne, deint
using Test

@testset "Deformula.jl" begin
    result = _deint(deformulaZeroToInf) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test result.s ≈ 1.0
    result = _deint(deformulaMinusOneToOne) do x
        y = x + 1
        0.5*y^-0.5*exp(-y^0.5)
    end
    @test result.s ≈ 0.7568830631028768
end

@testset "Deformula2" begin
    result = _deint(deformulaZeroToInf) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test sum(result.w) * result.h ≈ 1.0
end

@testset "Deformula3" begin
    result = deint(0.0, Inf64) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test sum(result.w) * result.h ≈ 1.0
end

@testset "Deformula4" begin
    result = deint(0.0, 10.0) do x
        if 0 <= x <= 10
            1
        else
            0
        end
    end
    println(result.x)
    @test result.s ≈ 10.0
end

@testset "Deformula5" begin
    result = deint(1.0, Inf64) do x
        1.0/sqrt(2.0*pi) * exp(-(x - 1.0)^2/2.0)
    end
    println(result.x)
    @test result.s ≈ 0.5
end

@testset "Deformula6" begin
    dfm(x; m, lambda) = lambda * x^(m-1)* exp(-lambda * x^m) + lambda * m * x^(m-1) * log(x) * exp(-lambda * x^m) + lambda * m * x^(m-1) * exp(-lambda * x^m) * (-lambda * x^m * log(x))
    result = deint(0.0, Inf64) do x
        dfm(x, m=2.0, lambda=1.0)
    end
    println(length(result.x))
    println(result.s)
    # @test result.s ≈ 0.5
end
