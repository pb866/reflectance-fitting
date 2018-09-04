# Try some tests of the Optim package to learn more about it.

using Optim
import PyPlot
const plt = PyPlot
using Printf
using Test
using ALS
import Random

# example from documentation
function rosenbrock()
    @testset "Rosenbrock with defaults from 0.0" begin
    f(x) = (1.0-x[1])^2+100.0*(x[2]-x[1])^2
    x0 = [0.0,0.0]
    res = optimize(f,x0)
    par = Optim.minimizer(res)
    @test isapprox(par[1], 1.0, atol=2e-5)
    @test isapprox(par[2], 1.0, atol=2.e-5)
    end
end

function rosenbrock_at_min()
    @testset "Rosenbrock with defaults from 1.0" begin
    f(x) = (1.0-x[1])^2+100.0*(x[2]-x[1])^2
    x0 = [1.0,1.0]
    res = optimize(f,x0)
    par = Optim.minimizer(res)
    @test isapprox(par[1], 1.0, atol=2e-5)
    @test isapprox(par[2], 1.0, atol=2.e-5)
    end
end

function rBFGS()
    @testset "Rosenbrock with BFGS from 0.0" begin
    f(x) = (1.0-x[1])^2+100.0*(x[2]-x[1])^2
    x0 = [0.0,0.0]
    res = optimize(f,x0,BFGS(), autodiff = :forward)
    par = Optim.minimizer(res)
    @test Optim.converged(res)
    @test Optim.iterations(res) > 0
    # @test Optim.x_converged(res)
    # @test Optim.f_converged(res)
    @test Optim.g_converged(res)
    @test isapprox(par[1], 1.0, atol=2e-5)
    @test isapprox(par[2], 1.0, atol=2.e-5)
    end
end

function rBFGS_at_min()
    @testset "Rosenbrock with with BFGS from 1.0" begin
    f(x) = (1.0-x[1])^2+100.0*(x[2]-x[1])^2
    x0 = [1.0,1.0]
    res = optimize(f,x0,BFGS())
    par = Optim.minimizer(res)
    @test Optim.converged(res)
    # @test Optim.iterations(res) > 0
    # @test Optim.x_converged(res)
    # @test Optim.f_converged(res)
    @test Optim.g_converged(res)
    @test isapprox(par[1], 1.0, atol=2e-5)
    @test isapprox(par[2], 1.0, atol=2.e-5)
    end
end

function refl()
    @testset "Nelder Mead exact" begin
        lam = 15.0
        AlIndex = Index("Al")
        alndx = AlIndex[lam]
        SiO2Index = Index("SiO2")
        sio2ndx = SiO2Index[lam]
        n = [1.0, alndx, sio2ndx]
        t = [0.0, 20.0, 0.0]
        thetad = 20.0
        p = parratt(n, t, thetad, lam)
        @test isapprox(p, 0.01260649028, atol = 1e-5)
        theta = 1.0:2.0:80.0
        rfl = [parratt(n, t, th, lam) for th in theta]
        function f(x,b)
            n=[1.0, b[1]+b[2]*im, sio2ndx]
            [parratt(n, t, th, lam) for th in x]
        end
        ssq(b) = sum(abs2,rfl-f(theta,b))
        x0 = [real(alndx), imag(alndx)]
        res = optimize(ssq, x0,
            Optim.Options(x_tol = 1e-10,
                f_tol = 1e-10,
                iterations = 10_000))
        par = Optim.minimizer(res)
        @test Optim.converged(res)
        @test isapprox(par[1], x0[1], atol=1e-7)
        @test isapprox(par[2], x0[2], atol=1e-7)
        @test Optim.minimum(res) < 1e-7
    end
end

function refl_BFGS()
    @testset "BFGS exact" begin
        lam = 15.0
        AlIndex = Index("Al")
        alndx = AlIndex[lam]
        SiO2Index = Index("SiO2")
        sio2ndx = SiO2Index[lam]
        n = [1.0, alndx, sio2ndx]
        t = [0.0, 20.0, 0.0]
        thetad = 20.0
        p = parratt(n, t, thetad, lam)
        @test isapprox(p, 0.01260649028, atol = 1e-5)
        theta = 1.0:2.0:80.0
        rfl = [parratt(n, t, th, lam) for th in theta]
        function f(x,b)
            n=[1.0, b[1]+b[2]*im, sio2ndx]
            [parratt(n, t, th, lam) for th in x]
        end
        ssq(b) = sum(abs2,rfl-f(theta,b))
        x0 = [real(alndx), imag(alndx)]
        res = optimize(ssq, x0, BFGS(),
            Optim.Options(x_tol = 1e-10,
                f_tol = 1e-10,
                iterations = 10_000))
        par = Optim.minimizer(res)
        @test Optim.converged(res)
        @test isapprox(par[1], x0[1], atol=1e-7)
        @test isapprox(par[2], x0[2], atol=1e-7)
        @test Optim.minimum(res) < 1e-7
    end
end

function refl_NM_noise()
    @testset "Nelder Mead noise" begin
        lam = 15.0
        AlIndex = Index("Al")
        alndx = AlIndex[lam]
        SiO2Index = Index("SiO2")
        sio2ndx = SiO2Index[lam]
        n = [1.0, alndx, sio2ndx]
        t = [0.0, 20.0, 0.0]
        thetad = 20.0
        p = parratt(n, t, thetad, lam)
        @test isapprox(p, 0.01260649028, atol = 1e-5)
        theta = 1.0:2.0:80.0
        rfl = [parratt(n, t, th, lam) for th in theta]
        Random.seed!(199382721)
        rfl .*= (1.0 .+ randn(length(theta))/10)
        function f(x,b)
            n=[1.0, b[1]+b[2]*im, sio2ndx]
            [parratt(n, t, th, lam) for th in x]
        end
        ssq(b) = sum(abs2,rfl-f(theta,b))
        x0 = [real(alndx), imag(alndx)]
        res = optimize(ssq, x0,
            Optim.Options(x_tol = 1e-10,
                f_tol = 1e-10,
                iterations = 10_000))
        par = Optim.minimizer(res)
        @test Optim.converged(res)
        @test isapprox(par[1], x0[1], rtol=0.01)
        @test isapprox(par[2], x0[2], rtol=0.08)
        @test Optim.minimum(res) < 0.02
    end
end

function refl_BFGS_noise()
    @testset "BFGS with noise" begin
        lam = 15.0
        AlIndex = Index("Al")
        alndx = AlIndex[lam]
        SiO2Index = Index("SiO2")
        sio2ndx = SiO2Index[lam]
        n = [1.0, alndx, sio2ndx]
        t = [0.0, 20.0, 0.0]
        thetad = 20.0
        p = parratt(n, t, thetad, lam)
        @test isapprox(p, 0.01260649028, atol = 1e-5)
        theta = 1.0:2.0:80.0
        rfl = [parratt(n, t, th, lam) for th in theta]
        Random.seed!(199382721)
        rfl .*= (1.0 .+ randn(length(theta))/10)
        function f(x,b)
            n=[1.0, b[1]+b[2]*im, sio2ndx]
            [parratt(n, t, th, lam) for th in x]
        end
        ssq(b) = sum(abs2,rfl-f(theta,b))
        x0 = [real(alndx), imag(alndx)]
        res = optimize(ssq, x0, BFGS(),
            Optim.Options(x_tol = 1e-10,
                f_tol = 1e-10,
                iterations = 10_000))
        par = Optim.minimizer(res)
        @test Optim.converged(res)
        @test isapprox(par[1], x0[1], rtol=0.01)
        @test isapprox(par[2], x0[2], rtol=0.08)
        @test Optim.minimum(res) < 0.02
    end
end

function refl_inexact()
    @testset "Inexact start" begin
        lam = 15.0
        AlIndex = Index("Al")
        alndx = AlIndex[lam]
        SiO2Index = Index("SiO2")
        sio2ndx = SiO2Index[lam]
        n = [1.0, alndx, sio2ndx]
        t = [0.0, 20.0, 0.0]
        thetad = 20.0
        p = parratt(n, t, thetad, lam)
        @test isapprox(p, 0.01260649028, atol = 1e-5)
        theta = 1.0:2.0:80.0
        rfl = [parratt(n, t, th, lam) for th in theta]
        Random.seed!(199382721)
        rfl .*= (1.0 .+ randn(length(theta))/10)
        function f(x,b)
            n=[1.0, b[1]+b[2]*im, sio2ndx]
            t = [0.0, b[3], 0.0]
            [parratt(n, t, th, lam) for th in x]
        end
        ssq(b) = sum(abs2,rfl-f(theta,b))
        x0 = [real(alndx)*0.8, imag(alndx)*1.2, 30.0]
        res = optimize(ssq, x0, BFGS(),
            Optim.Options(x_tol = 1e-10,
                f_tol = 1e-10,
                iterations = 10_000))
        par = Optim.minimizer(res)
        @test Optim.converged(res)
        @test isapprox(par[1], real(alndx), rtol=0.01)
        @test isapprox(par[2], imag(alndx), rtol=0.08)
        @test isapprox(par[3], 20.0, atol=0.2)
        @test Optim.minimum(res) < 0.02
        # Plot the result
        xx = 1.0:0.1:80.0
        yy = f(xx,par)
        plt.figure()
        plt.semilogy(theta, rfl, ".", xx,yy,"-")
    end
end

function dotests()
    @testset "Rosenbrock" begin
        rosenbrock()
        rosenbrock_at_min()
        rBFGS()
        rBFGS_at_min();
    end
    @testset "refl" begin
        refl();
        refl_BFGS();
        refl_NM_noise()
        refl_BFGS_noise()
        refl_inexact()
    end
end

dotests();
