# uncert.jl
# Code to test uncertainties.tex

# Simple linear fit
using LinearAlgebra
using Printf

function lfit()
    npts = 10000
    xpts = [(i-1.0)/(npts-1) for i=1:npts]
    m = 1.5
    b = 3.0
    ypts = b .+ m.*xpts .+ randn(npts)*0.2
    A = Matrix{Union{Missing, Float64}}(missing, 2, 2)
    y = Array{Union{Missing, Float64}}(missing, 2)
    A[1,1] = npts
    A[1,2] = sum(xpts)
    A[2,1] = A[1,2]
    A[2,2] = sum(xpts.*xpts)
    y[1] = sum(ypts)
    y[2] = sum(ypts.*xpts)
    b = A\y
    J = hcat(ones(npts),xpts)
    yfit = b[1] .+ xpts.*b[2]
    res = ypts .- yfit
    sigmay = sqrt.(sum(abs2,res)/(npts-2))
    cov = inv(J'*J)
    sigmai = sigmay.*sqrt.(diag(cov))
    [b sigmai]
end

# repeat 1000 times for averages
let
    psum = zeros(2,2)
    sumsq = zeros(2,2)
    trials = 1000
    for i=1:trials
        ft = lfit()
        psum += ft
        sumsq += ft.*ft
    end
    pbar = psum./trials
    pbarsq = sumsq./trials
    sigma = sqrt.(pbarsq[:,1] .- pbar[:,1].^2)
    @printf("average intercept = %.3f +/- %.4f\n",
        pbar[1,1], pbar[1,2])
    @printf("average slope = %.3f +/- %.4f\n",
        pbar[2,1], pbar[2,2])
    @printf("computed uncertainty in intercept: %.2e\n",
        sigma[1])
    @printf("computed uncertainty in slope: %.2e\n",
        sigma[2])
end
