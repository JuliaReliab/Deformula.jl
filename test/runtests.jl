using Deformula
using Test

@testset "Deformula.jl" begin
    result = deint(deformulaZeroToInf) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test result[1] ≈ 1.0
    result = deint(deformulaMinusOneToOne) do x
        y = x + 1
        0.5*y^-0.5*exp(-y^0.5)
    end
    @test result[1] ≈ 0.7568830631028768
end

@testset "Deformula2" begin
    result = deint(deformulaZeroToInf) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test sum(result[4]) ≈ 1.0
end
