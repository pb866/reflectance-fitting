# uncert.jl
# Code to test uncertainties.tex

# Nonlinear Fit
using LsqFit
using Printf
using LinearAlgebra

function nlfit()
    npts = 10000
    f(x,p) = p[1]./(p[2] .+ cos.(2*pi .* x ./p[3]))
    xpts = [(i-1.0)/(npts-1) for i=1:npts]
    p0 = [0.5,1.35,0.3]
    ypts = f(xpts,p0)+randn(npts)/0.2
    wt = ones(npts)
    cf = curve_fit(f, xpts, ypts, wt, p0)
    sigmay = sqrt.(sum(abs2,cf.resid)/(cf.dof))
    J = cf.jacobian
    cov = inv(J'*J)
    sigmai = sigmay.*sqrt.(diag(cov))
    [cf.param sigmai]
end

# repeat 1000 times for averages
let
    psum = zeros(3,2)
    sumsq = zeros(3,2)
    trials = 1000
    for i=1:trials
        ft = nlfit()
        psum += ft
        sumsq += ft.*ft
    end
    pbar = psum./trials
    pbarsq = sumsq./trials
    sigma = sqrt.(pbarsq[:,1] .- pbar[:,1].^2)
    @printf("average b1 = %.3f +/- %.4f\n",
        pbar[1,1], pbar[1,2])
    @printf("average b2 = %.3f +/- %.4f\n",
        pbar[2,1], pbar[2,2])
    @printf("average b3 = %.3f +/- %.4f\n",
            pbar[3,1], pbar[3,2])
    @printf("computed uncertainty in b1: %.2e\n",
        sigma[1])
    @printf("computed uncertainty in b2: %.2e\n",
        sigma[2])
        @printf("computed uncertainty in b3: %.2e\n",
            sigma[3])
end
