"""
    LsqFit Module

Steve Turley's implementation of the LsqFit Julia library module.
Defines and exports the following functions:
* curve_fit
* linear_fit
* poly_fit
* standard_error

```curve_fit``` returns its results in a ```CurveFitResult``` structure.
"""
module LsqFit
# LsqFit replacement

export curve_fit, standard_error, linear_fit, poly_fit

using Optim
using Calculus
using LinearAlgebra
using Test
import Random
import LineSearches

# if true, module does unit testing
const dotests = false

"""
    CurveFitResult(dof, param, resid, jacobian, converged, wt, mse)

concrete type holding the results from a call to ```curve_fit```
"""
struct CurveFitResult
    dof::Int
    param::Vector{Float64}
    resid::Vector{Float64}
    jacobian::Matrix{Float64}
    converged::Bool
    wt::Array{Float64}
    mse::Float64
end

"""
    curve_fit(model, xpts, ypts, wt, p0[, optimizer])

# Parameters
* ``model``: function with two arguments, x (array of independent data) and p (parameters to be fit). Should return a vector with the function value at the x points.
* ``xpts``: independent variable for the data
* ``ypts``: dependent variable with the data
* ``wt``: array of weights for the fit. Function fit is sum(abs2, wt.*(model(x,p)-ypts))
* ``p0``: initial guess for fit parameters
* ``optimizer``: optional command to specify the optimizer to use. Choices are :BFGS and :NelderMead.
"""
function curve_fit(model::Function, xpts::AbstractArray, ypts::AbstractArray,
        wt::AbstractArray, p0::Vector; optimizer::AbstractString="NelderMead",
        linesearch=LineSearches.HagerZhang())
    # construct a weighted cost function, with a vector weight for each ydata
    # for example, this might be wt = 1/sigma where sigma is some error term
    f(p) = wt .* ( model(xpts, p) - ypts )
    ssq(p) = sum(abs2, f(p))
    options = Optim.Options(x_tol = 1e-10, f_tol = 1e-10, iterations = 10_000)
    if optimizer=="BFGS"
        results = optimize(ssq, p0, BFGS(linesearch=linesearch), options)
    elseif optimizer=="NelderMead"
        results = optimize(ssq, p0, NelderMead(linesearch=linesearch), options)
    else
        error("specified unknown optimizer")
    end
    # Fill the result with useful data
    dof = length(xpts) - length(p0)
    p = Optim.minimizer(results)
    res = f(p)
    g = Calculus.jacobian(f)
    jac = g(p)
    conv = Optim.converged(results)
    mse = sum(abs2, res)/dof
    CurveFitResult(dof, p, res, jac, conv, wt, mse)
end

"""
    linear_fit(model, xpts, ypts)

# Parameters
* basis: array of basis function to use in the fit. They should be functions of a vector of points x which return a vector of points y
* xpts: independent variable for the data
* ypts: dependent variable with the data
"""
function linear_fit(basis::AbstractArray, xpts::AbstractArray,
        ypts::AbstractArray)
    order=length(basis)
    A=zeros(order,order)
    y=zeros(order)
    for i=1:order
        for j=1:order
            A[i,j] = sum(basis[i](xpts).*basis[j](xpts))
        end
        y[i] = sum(ypts.*basis[i](xpts))
    end
    b = A\y
    npts = length(xpts)
    dof = npts-order
    yf = zeros(npts)
    for i=1:order
        yf .+= b[i].*basis[i](xpts)
    end
    res = yf .- ypts
    jac=Matrix{Union{Missing, Float64}}(missing,order,npts)
    for i=1:order
        jac[i,:] = basis[i](xpts)
    end
    conv = true
    wt = ones(npts)
    mse = sum(abs2, res)/dof
    CurveFitResult(dof, b, res, jac', conv, wt, mse)
end

"""
    poly_fit(order, xpts, ypts)

Fit a polynomial to the data in xpts and ypts.
# Parameters
* order: highest order monomial to use in the fit
* xpts: independent variable for the data
* ypts: dependent variable with the data
"""
function poly_fit(order::Int, xpts::AbstractArray, ypts::AbstractArray)
    lf = [x->x.^i for i=0:order]
    linear_fit(lf, xpts, ypts)
end

function estimate_covar(fit::CurveFitResult)
    # computes covariance matrix of fit parameters
    J = fit.jacobian
    inv(J'*J)
end

"""
    standard_error(fit, [rtol, atol])

Estimate the standard deviation in the fit parameters returned by
```curve_fit```.
# Parameters
* fit: CurveFitResult returned by ```curve_fit```
* rtol: relative error (ignored in this implementation)
* atol: absolute error (ignored in this implementation)
"""
function standard_error(fit::CurveFitResult; rtol::Real=NaN, atol::Real=0)
    # computes standard error of estimates from
    #   fit   : a CurveFitResult from a curve_fit()
    covar = estimate_covar(fit)
    # then the standard errors are given by the sqrt of the diagonal
    vars = diag(covar)
    # Take the absolute value to be safe
    sqrt.(abs.(vars*fit.mse))
end

function test_one()
    @testset "Single Fit" begin
    model(x,p) = exp.(-p[1].*x).*(1.1.+0.8.*cos.(p[2].*x))/2
    p=[1,pi];
    x=range(0,stop=10,length=200);
    y=model(x,p)
    fval = 0.2
    ypts=model(x,p).*(1.0.+fval.*randn(length(x)));

    w=1.0./(fval.*y);
    mfit = curve_fit(model, x, ypts, w, p);

    @test mfit.dof == 198
    @test isapprox(mfit.param[1], 1.0, atol=1e-2)
    @test isapprox(mfit.param[2], pi, atol=1e-2)
    sigma = standard_error(mfit)
    @test isapprox(sigma[1], 2.44e-3, atol = 2e-4)
    @test isapprox(sigma[2], 3.64e-3, atol = 2e-4)
    @test mfit.converged # doesn't always pass
    @test isapprox(mfit.mse, 1.0, atol = 0.25)
    end;
end

# Checking uncertainties
# Do multiple fits and compare estimated uncertainties with
# actual variations in the data.
function afit(trials::Int)
    model(x,p) = exp.(-p[1].*x).*(1.1.+0.8.*cos.(p[2].*x))/2
    p=[1,pi];
    x=range(0,stop=10,length=200);
    y=model(x,p)
    fval = 0.2
    psum=[0.0,0.0]
    p2sum=[0.0,0.0]
    vsum=[0.0,0.0]
    fvar = 0.2;
    for i = 1:trials
        ypts=model(x,p).*(1.0.+fvar.*randn(length(x)));
        w=1.0./(fvar.*y);
        fit = curve_fit(model, x, ypts, w, p);
        sigma = standard_error(fit)
        psum = psum + fit.param
        p2sum = p2sum + fit.param.^2
        vsum = vsum + sigma
    end
    pbar = psum/trials
    p2bar = p2sum/trials
    vbar = vsum/trials
    sigmap = sqrt.(p2bar-pbar.^2)
    (vbar, sigmap)
end

# Checking uncertainties
# Do multiple fits and compare estimated uncertainties with
# actual variations in the data.
function lfit(trials::Int)
    fl = [x->ones(length(x)), x->x, x->x.^2]
    Random.seed!(199382721)
    xpts = 1:10;
    lf = linear_fit(fl, xpts, ypts)
    b = lf.param
    std = standard_error(lf)
    psum=zeros(3)
    p2sum=zeros(3)
    vsum=zeros(3)
    fvar = 10.0;
    for i = 1:trials
        xnoise = randn(10)./fvar;
        ypts = xnoise .+ 5.0 .+ 3.0.*xpts .+ 0.1.*xpts.^2
        fit = linear_fit(fl, xpts, ypts);
        sigma = standard_error(fit)
        psum = psum + fit.param
        p2sum = p2sum + fit.param.^2
        vsum = vsum + sigma
    end
    pbar = psum/trials
    p2bar = p2sum/trials
    vbar = vsum/trials
    sigmap = sqrt.(p2bar-pbar.^2)
    (vbar, sigmap)
end

function test_linear()
    @testset "linear_fit check" begin
        fl = [x->ones(length(x)), x->x, x->x.^2]
        Random.seed!(199382721)
        xnoise = randn(10)./10.0;
        xpts = 1:10;
        ypts = xnoise .+ 5.0 .+ 3.0.*xpts .+ 0.1.*xpts.^2
        lf = linear_fit(fl, xpts, ypts)
        b = lf.param
        @test length(b) == 3
        @test abs(b[1]-5.0) < 0.03
        @test abs(b[2]-3.0) < 0.02
        @test abs(b[3]-0.1) < 0.002
        @test lf.mse < 0.01
        @test lf.mse > 0.009
        vb, sb = afit(1000)
        @test isapprox(vb[1], sb[1], rtol = 0.04)
        @test isapprox(vb[2], sb[2], rtol = 0.02)
    end;
end

function test_poly()
    @testset "poly_fit check" begin
        Random.seed!(199382721)
        xnoise = randn(10)./10.0;
        xpts = 1:10
        ypts = xnoise .+ 5.0 .+ 3.0.*xpts .+ 0.1.*xpts.^2
        lf = poly_fit(2, xpts, ypts)
        b = lf.param
        @test length(b) == 3
        @test abs(b[1]-5.0) < 0.03
        @test abs(b[2]-3.0) < 0.02
        @test abs(b[3]-0.1) < 0.002
        @test lf.mse < 0.01
        @test lf.mse > 0.009
        vb, sb = afit(1000)
        @test isapprox(vb[1], sb[1], rtol = 0.04)
        @test isapprox(vb[2], sb[2], rtol = 0.02)
    end;
end

if dotests
    @testset "LsqFit" begin
        test_one()
        @testset "Covariance Check" begin
            vb, sb = afit(1000)
            @test isapprox(vb[1], sb[1], rtol = 8e-2)
            @test isapprox(vb[2], sb[2], rtol = 0.05)
        end;
        test_linear()
        test_poly()
        end;
end

# Since the values were so similar, this must be a good way to do it.
# Note that I tried this with different sigmas (fvar) and they still
# agreed.
end
