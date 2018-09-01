# Simulate a t-test
# Start with one variable distributed about 0
using Printf
using Distributions

const doPlot = false

if doPlot
    import PyPlot
    const pp = PyPlot
    td = TDist(1)
    x=-3:0.05:3;
    tpdf = [pdf(td,y) for y in x];
    tcdf = [cdf(td,y) for y in x];

    pp.figure();
    pp.plot(x,tpdf);
    pp.title("T Distribution PDF")

    pp.figure()
    pp.plot(x,tcdf);
    pp.title("T Distribution CDF")
end

# Compare with Excel
smean = 4.694206958;
serr = 0.430462388;
tStat = smean/serr;
@printf("tStat = %f\n", tStat)
for df = 1:5
    td = TDist(df)
    tpdf = pdf(td, tStat)
    tcdf = cdf(td, tStat)
    @printf("For dof = %d pdf = %.5f, cdf = %.5f, 2(1-cdf) = %.5f\n",
        df, tpdf, tcdf, 2*(1-tcdf))
end

# Conclusion p value is 2(1-cdf) with npts-params
# degrees of freedom


# Recreate excel fits
noise=[0.712705059,
    1.037487891,
    0.292159257,
    -0.71952627,
    -1.187286216]
xpts = 1:5
m = 15
b = 3
ypts = noise .+ m.*xpts .+ b
using LsqFit
pf = poly_fit(1,xpts,ypts)
sigma = standard_error(pf)
@printf("fit intercept = %.5f +/- %.5f\n",
    pf.param[1], sigma[1])
@printf("fit slope = %.5f +/- %.5f\n",
        pf.param[2], sigma[2])

# Check t-Stat and P-value
tstat = pf.param./sigma
@printf("fit tStats = %.2f, %.1f\n", tstat...)
df = 3
td = TDist(df)
tcdf = [cdf(td, tstat[i]) for i=1:2]
pval = 2 .* (1 .-tcdf)
@printf("P-value = %.2e %.2e\n", pval...)


# uncertainties agree, check distribution

function tfit(n::Int)
    xpts = 1:5
    m = 15
    b = 3
    psum = zeros(2)
    psum2 = zeros(2)
    ssum = zeros(2)
    for i=1:n
        noise = randn(5);
        ypts = noise .+ m.*xpts .+ b
        pf = poly_fit(1, xpts, ypts)
        sigma = standard_error(pf)
        psum += pf.param
        psum2 += pf.param.^2
        ssum += sigma
    end
    pbar = psum./n
    pbar2 = psum2./n
    sbar = ssum./n
    pstd = sqrt.(pbar2-pbar.^2)
    @printf("Average intercept = %.5f +/- %.5f\n",
        pbar[1], pstd[1])
    @printf("Average slope = %.5f +/- %.5f\n",
        pbar[2], pstd[2])
    println("Average predicted standard deviation = ",
        pstd)
end

tfit(10000)
