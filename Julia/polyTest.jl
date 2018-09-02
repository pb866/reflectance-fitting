# Testing polynomial fitting

using LsqFit
using LinearAlgebra
using PyPlot
const plt = PyPlot
using Printf

function printmat(str, a)
    println(str)
    for i=1:3
        @printf("%7.3f %7.3f %7.3f\n",a[i,:]...)
    end
    println()
end

function printemat(str, a)
    println(str)
    for i=1:3
        @printf("%10.2e %10.2e %10.2e\n",a[i,:]...)
    end
    println()
end

npts = 1000
xpts = [(i-1)/(npts-1) for i=1:1000]
b0 = 2.3
b1 = 1.7
b2 = 1.4
ypts = b0 .+ b1.*xpts .+ b2.*xpts.^2+randn(npts)*0.2
pf = poly_fit(2,xpts,ypts)
yfit = pf.param[1] .+ pf.param[2].*xpts .+
    pf.param[3].*xpts.^2
plt.figure()
plt.plot(xpts, ypts, ".", xpts, yfit, "-")
plt.title("Data and Fit for Polynomial Fit")
J = pf.jacobian
cov = pf.mse.*inv(J'*J)
println("polynomial\n")
@printf("b1 = %.3f b2 = %.3f b3 = %.3f\n\n", pf.param...)
printemat("cov:",cov)
ev = LinearAlgebra.eigen(cov)
printmat("eigenvectors:",ev.vectors)

println("\northgonal polynomials\n")
basis = [x->ones(length(x)), x->2 .* x .- 1, x->6 .* x.^2 .- 6 .* x .+ 1]
lf = linear_fit(basis, xpts, ypts)
@printf("b1 = %.3f b2 = %.3f b3 = %.3f\n\n", lf.param...)
J = lf.jacobian
cov = lf.mse.*inv(J'*J)
printemat("cov:",cov)
ev = LinearAlgebra.eigen(cov)
printmat("eigenvectors:",ev.vectors)
