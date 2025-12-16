import DEQuadrature: _deint, deformulaZeroToInf, deformulaMinusOneToOne, deint
using Test

@testset "DEQuadrature.jl" begin
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

# Additional tests

@testset "FiniteIntervalGeneral" begin
    # Integrate constant 1 over a non-symmetric finite interval
    a, b = -2.5, 3.7
    result = deint(a, b) do x
        1.0
    end
    @test result.s ≈ (b - a) atol=1e-8 rtol=1e-8
    @test issorted(result.t)
end

@testset "BigFloatSemiInfinite" begin
    setprecision(256) do
        result = deint(BigFloat(0), BigFloat(Inf); reltol=BigFloat(1e-20)) do x
            exp(-x)
        end
        @test result.s ≈ BigFloat(1) rtol=BigFloat(1e-18)
    end
end

@testset "BigFloatFiniteInterval" begin
    setprecision(256) do
        a = BigFloat(-2)
        b = BigFloat(3)
        result = deint(a, b; reltol=BigFloat(1e-20)) do x
            one(BigFloat)
        end
        @test result.s ≈ (b - a) rtol=BigFloat(1e-12)
        @test issorted(result.t)
    end
end

@testset "MaxIterWarning" begin
    # Force a warning by setting maxiter very small; ensure it doesn't throw
    @test_logs (:warn,) begin
        result = _deint(deformulaZeroToInf; maxiter=1) do x
            exp(-x)
        end
        @test result.s > 0
    end
end

@testset "SemiInfiniteShift" begin
    # Check [lower, Inf) with a shift works; mean=2 normal -> P(X>=2)=0.5
    result = deint(2.0, Inf64) do x
        1.0 / sqrt(2.0*pi) * exp(-((x - 2.0)^2) / 2.0)
    end
    @test result.s ≈ 0.5 atol=1e-8 rtol=1e-8
end

@testset "Deformula6" begin
    dfm(x; m, lambda) = lambda * x^(m-1)* exp(-lambda * x^m) + lambda * m * x^(m-1) * log(x) * exp(-lambda * x^m) + lambda * m * x^(m-1) * exp(-lambda * x^m) * (-lambda * x^m * log(x))
    # This integrand is the partial derivative w.r.t. m of the Weibull PDF λ m x^{m-1} e^{-λ x^m}.
    # Therefore, ∫_0^∞ dfm(x; m, λ) dx = ∂/∂m ∫_0^∞ PDF dx = ∂/∂m 1 = 0.
    # We assert the integral is approximately 0 within absolute tolerance.
    result = deint(0.0, Inf64; reltol=1e-8, abstol=1e-12, d=8, maxiter=20) do x
        dfm(x, m=2.0, lambda=1.0)
    end
    @test isfinite(result.s)
    @test result.s ≈ 0.0 atol=1e-8
end
