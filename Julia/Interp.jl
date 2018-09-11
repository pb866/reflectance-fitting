module Interp
# Code for interpolation for various orders
using LinearAlgebra
using Test
import Base.length
# import Base.getindex

export CubicSpline, interp, slope, slope2, pchip, pchip2, pchip3

const unitTests = false
const graphicsTests = false
const bumpTests = false
if graphicsTests | bumpTests
    import PyPlot # needed if graphicsTests is true
    const plt = PyPlot
end
const eps = 1e-4 # rel error allowed on extrapolation

"""
    CubicSpline(x,a,b,c,d)

concrete type for holding the data needed
    to do a cubic spline interpolation
"""
struct CubicSpline
    x::Union{Array{Float64,1},
        StepRangeLen{Float64,
            Base.TwicePrecision{Float64},
            Base.TwicePrecision{Float64}}}
    a::Array{Float64,1}
    b::Array{Float64,1}
    c::Array{Float64,1}
    d::Array{Float64,1}
end

"""
    PCHIP(x,a,b,c,d)

concrete type for holding the data needed
    to do a piecewise continuous hermite interpolation
"""
struct PCHIP
    x::Union{Array{Float64,1},
        StepRangeLen{Float64,
            Base.TwicePrecision{Float64},
            Base.TwicePrecision{Float64}}}
    y::Array{Float64,1}
    d::Array{Float64,1}
    h::Array{Float64,1}
end

"""
    CubicSpline(x,y)

Creates the CubicSpline structure needed for cubic spline
interpolation

# Arguments
- `x`: an array of x values at which the function is known
- `y`: an array of y values corresonding to these x values
"""
function CubicSpline(x::Array{Float64,1}, y::Array{Float64,1})
    len = size(x,1)
    if len<3
        error("CubicSpline requires at least three points for interpolation")
    end
    # Pre-allocate and fill columns and diagonals
    yy = zeros(len)
    dl = zeros(len-1)
    du = zeros(len-1)
    dd = zeros(len)
    alpha = x[2:len].-x[1:len-1]
    yy[1] = 3*(y[2]-y[1])/alpha[1]^2
    du[1] = 1/alpha[1]
    dd[1] = 2/alpha[1]
    for i=2:len-1
        yy[i] = 3*(y[i+1]/alpha[i]^2+
            y[i]*(alpha[i-1]^(-2)-alpha[i]^(-2))-
            y[i-1]/alpha[i-1]^2)
        dl[i-1] = 1/alpha[i-1]
        du[i] = 1/alpha[i]
        dd[i] = 2*(1/alpha[i-1]+1/alpha[i])
    end
    yy[len] = 3*(y[len]-y[len-1])/alpha[len-1]^2
    dl[len-1] = 1/alpha[len-1]
    dd[len] = 2/alpha[len-1]
    # Solve the tridiagonal system for the derivatives D
    dm = Tridiagonal(dl,dd,du)
    D = dm\yy
    # fill the arrays of spline coefficients
    a = y[1:len-1]  # silly but makes the code more transparent
    b = D[1:len-1]  # ditto
    c = 3 .*(y[2:len].-y[1:len-1])./alpha[1:len-1].^2 .-
        2*D[1:len-1]./alpha[1:len-1].-D[2:len]./alpha[1:len-1]
    d = 2 .*(y[1:len-1].-y[2:len])./alpha[1:len-1].^3 .+
        D[1:len-1]./alpha[1:len-1].^2 .+
        D[2:len]./alpha[1:len-1].^2
    CubicSpline(x,a,b,c,d)
end

function CubicSpline(x::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, y::Array{Float64,1})
    len = length(x)
    if len<3
        error("CubicSpline requires at least three points for interpolation")
    end
    # Pre-allocate and fill columns and diagonals
    yy = zeros(len)
    dl = ones(len-1)
    dd = 4.0 .* ones(len)
    dd[1] = 2.0
    dd[len] = 2.0
    yy[1] = 3*(y[2]-y[1])
    for i=2:len-1
        yy[i] = 3*(y[i+1] - y[i-1])
    end
    yy[len] = 3*(y[len]-y[len-1])
    # Solve the tridiagonal system for the derivatives D
    dm = Tridiagonal(dl,dd,dl)
    D = dm\yy
    # fill the arrays of spline coefficients
    a = y[1:len-1]  # silly but makes the code more transparent
    b = D[1:len-1]  # ditto
    c = 3 .*(y[2:len].-y[1:len-1]).-2*D[1:len-1].-D[2:len]
    d = 2 .*(y[1:len-1].-y[2:len]).+D[1:len-1].+D[2:len]
    CubicSpline(x,a,b,c,d)
end

# This version of pchip uses the mean value of the slopes
# between data points on either side of the interpolation point
"""
    pchip(x,y)

Creates the PCHIP structure needed for piecewise
    continuous cubic spline interpolation

# Arguments
- `x`: an array of x values at which the function is known
- `y`: an array of y values corresonding to these x values
"""
function pchip(x::Array{Float64,1}, y::Array{Float64,1})
    len = size(x,1)
    if len<3
        error("PCHIP requires at least three points for interpolation")
    end
    h = x[2:len].-x[1:len-1]
    # Pre-allocate and fill columns and diagonals
    d = zeros(len)
    d[1] = (y[2]-y[1])/h[1]
    for i=2:len-1
        d[i] = (y[i+1]/h[i]+y[i]*(1/h[i-1]-1/h[i])-y[i-1]/h[i-1])/2
    end
    d[len] = (y[len]-y[len-1])/h[len-1]
    PCHIP(x,y,d,h)
end

# PCHIP with quadratic fit to determine slopes
function pchip2(x::Array{Float64,1}, y::Array{Float64,1})
    len = size(x,1)
    if len<3
        error("PCHIP requires at least three points for interpolation")
    end
    h = x[2:len].-x[1:len-1]
    # Pre-allocate and fill columns and diagonals
    d = zeros(len)
    d[1] = (y[2]-y[1])/h[1]
    for i=2:len-1
        d[i] = (y[i]-y[i-1])*h[i]/(h[i-1]*(h[i-1]+h[i])) +
            (y[i+1]-y[i])*h[i-1]/(h[i]*(h[i-1]+h[i]))
    end
    d[len] = (y[len]-y[len-1])/h[len-1]
    PCHIP(x,y,d,h)
end

# Real PCHIP
function pchip3(x::Array{Float64,1}, y::Array{Float64,1})
    len = size(x,1)
    if len<3
        error("PCHIP requires at least three points for interpolation")
    end
    h = x[2:len].-x[1:len-1]
    del = (y[2:len].-y[1:len-1])./h
    # Pre-allocate and fill columns and diagonals
    d = zeros(len)
    d[1] = del[1]
    for i=2:len-1
        if del[i]*del[i-1] < 0
            d[i] = 0
        else
            d[i] = (del[i]+del[i-1])/2
        end
    end
    d[len] = del[len-1]
    for i=1:len-1
        if del[i] == 0
            d[i] = 0
            d[i+1] = 0
        else
            alpha = d[i]/del[i]
            beta = d[i+1]/del[i]
            if alpha^2+beta^2 > 9
                tau = 3/sqrt(alpha^2+beta^2)
                d[i] = tau*alpha*del[i]
                d[i+1] = tau*beta*del[i]
            end
        end
    end
    PCHIP(x,y,d,h)
end

"""
    interp(cs::CubicSpline, v::Float)

Interpolate to the value corresonding to v

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
cs = CubicSpline(x,y)
v = interp(cs, 1.2)
```
"""
function interp(cs::CubicSpline, v::Float64)
    # Find v in the array of x's
    if (v<cs.x[1]) | (v>cs.x[length(cs.x)])
        error("Extrapolation not allowed")
    end
    segment = region(cs.x, v)
    if typeof(cs.x)==Array{Float64,1}
        # irregularly spaced points
        t = v-cs.x[segment]
    else
        # regularly spaced points
        t = (v-cs.x[segment])/(cs.x[segment+1]-cs.x[segment])
    end
    cs.a[segment] + t*(cs.b[segment] + t*(cs.c[segment] + t*cs.d[segment]))
end

function interp(pc::PCHIP, v::Float64)

    if v*(1+eps)<first(pc.x)
        error("Extrapolation not allowed, $v<$(first(pc.x))")
    end
    if v*(1-eps)>last(pc.x)
        error("Extrapolation not allowed, $v>$(last(pc.x))")
    end
    i = region(pc.x, v)
    phi(t) = 3*t^2 - 2*t^3
    psi(t) = t^3 - t^2
    H1(x) = phi((pc.x[i+1]-v)/pc.h[i])
    H2(x) = phi((v-pc.x[i])/pc.h[i])
    H3(x) = -pc.h[i]*psi((pc.x[i+1]-v)/pc.h[i])
    H4(x) = pc.h[i]*psi((v-pc.x[i])/pc.h[i])
    pc.y[i]*H1(v) + pc.y[i+1]*H2(v) + pc.d[i]*H3(v) + pc.d[i+1]*H4(v)
end

"""
    slope(cs::CubicSpline, v::Float)

Derivative at the point corresonding to v

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
cs = CubicSpline(x,y)
v = slope(cs, 1.2)
```
"""
function slope(cs::CubicSpline, v::Float64)
    # Find v in the array of x's
    if (v<cs.x[1]) | (v>cs.x[length(cs.x)])
        error("Extrapolation not allowed")
    end
    segment = region(cs.x, v)
    if typeof(cs.x)==Array{Float64,1}
        # irregularly spaced points
        t = v-cs.x[segment]
    else
        # regularly spaced points
        t = (v-cs.x[segment])/(cs.x[segment+1]-cs.x[segment])
    end
    cs.b[segment] + t*(2*cs.c[segment] + t*3*cs.d[segment])
end

"""
    slope(pc::PCHIP, v::Float)

Derivative at the point corresonding to v

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
pc = pchip(x,y)
v = slope(pc, 1.2)
```
"""
function slope(pc::PCHIP, v::Float64)
    # Find v in the array of x's
    if (v<pc.x[1]) | (v>pc.x[length(pc.x)])
        error("Extrapolation not allowed")
    end
    i = region(pc.x, v)
    phip(t) = 6*t - 6*t^2
    psip(t) = 3*t^2 - 2*t
    H1p(x) = -phip((pc.x[i+1]-v)/pc.h[i])/pc.h[i]
    H2p(x) = phip((v-pc.x[i])/pc.h[i])/pc.h[i]
    H3p(x) = psip((pc.x[i+1]-v)/pc.h[i])
    H4p(x) = psip((v-pc.x[i])/pc.h[i])
    pc.y[i]*H1p(v) + pc.y[i+1]*H2p(v) + pc.d[i]*H3p(v) + pc.d[i+1]*H4p(v)
end

"""
    slope2(cs::CubicSpline, v::Float)

Second derivative at the point corresonding to v

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
cs = CubicSpline(x,y)
v = slope2(cs, 1.2)
```
"""
function slope2(cs::CubicSpline, v::Float64)
    # Find v in the array of x's
    if (v<cs.x[1]) | (v>cs.x[length(cs.x)])
        error("Extrapolation not allowed")
    end
    segment = region(cs.x, v)
    if typeof(cs.x)==Array{Float64,1}
        # irregularly spaced points
        t = v-cs.x[segment]
    else
        # regularly spaced points
        t = (v-cs.x[segment])/(cs.x[segment+1]-cs.x[segment])
    end
    2*cs.c[segment] + 6*t*cs.d[segment]
end

function region(x::Array{Float64,1}, v::Float64)
    # Binary search
    len = size(x,1)
    li = 1
    ui = len
    mi = div(li+ui,2)
    done = false
    while !done
        if v<x[mi]
            ui = mi
            mi = div(li+ui,2)
        elseif v>x[mi+1]
            li = mi
            mi = div(li+ui,2)
        else
            done = true
        end
        if mi == li
            done = true
        end
    end
    mi
end

function region(x::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, y::Float64)
    min(trunc(Int,(y-first(x))/step(x)),length(x)-2) + 1
end

function regular_tests()
    @testset "regular interpolation" begin
    # Test not enough points exception
    x = range(1.0, stop=2.0, length=2)
    y = [2.0, 4.0]
    @test_throws ErrorException CubicSpline(x,y)
    x = range(1.0, stop=3.25, length=4)
    y = [1.5, 3.0, 3.7, 2.5]
    cs = CubicSpline(x,y)
    @test_throws ErrorException interp(cs, 0.0)
    @test_throws ErrorException interp(cs, 4.0)
    # Check region
    @test region(x, 1.0) == 1
    @test region(x, 1.2) == 1
    @test region(x, 3.25) == 3
    @test region(x, 2.0) == 2
    @test region(x, 2.8) == 3
    # Check spline at knots
    @test interp(cs, 1.0) == 1.5
    @test interp(cs, 1.75) == 3.0
    @test isapprox(interp(cs, 3.25), 2.5, atol=1e-14)
    # Check spline with unit spacing of knots
    x = range(0.0, stop=4.0, length=5)
    y = sin.(x)
    cs = CubicSpline(x,y)
    dy = cos.(x)
    for i = 1:4
        @test cs.a[i] == y[i]
        @test isapprox(cs.a[i] + cs.b[i] + cs.c[i] + cs.d[i], y[i+1], atol=1.e-12)
        @test isapprox(cs.b[i], dy[i], atol=0.08)
        @test isapprox(cs.b[i] + 2*cs.c[i] + 3*cs.d[i], dy[i+1], atol=0.25)
    end
    end;
end

function irregular_tests()
    x = range(0.0, stop=4.0, length=5)
    y = sin.(x)
    cs = CubicSpline(x,y)
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = sin.(x)
    csi = CubicSpline(x,y)
    @testset "irregular interpolation" begin
    # Test not enough points exception
    x = [1.0, 2.0]
    y = [2.0, 4.0]
    @test_throws ErrorException CubicSpline(x,y)
    x = [0.2, 1.4, 3.8, 5.7]
    y = [1.5, 3.0, 3.7, 2.5]
    csi = CubicSpline(x,y)
    @test_throws ErrorException interp(csi, 0.0)
    @test_throws ErrorException interp(csi, 6.0)
    # Check region
    @test region(x, 0.3) == 1
    @test region(x, 0.2) == 1
    @test region(x, 5.7) == 3
    @test region(x, 2.1) == 2
    @test region(x, 4.0) == 3
    # Check spline at knots
    @test interp(csi, 0.2) == 1.5
    @test interp(csi, 1.4) == 3.0
    @test isapprox(interp(csi, 5.7), 2.5, atol=1e-14)
    # Check spline with unit spacing of knots
    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = sin.(x)
    csi = CubicSpline(x,y)
    for i = 1:4
        @test csi.a[i] == cs.a[i]
        @test csi.b[i] == cs.b[i]
        @test csi.c[i] == cs.c[i]
        @test csi.d[i] == cs.d[i]
        @test csi.a[i] == y[i]
        @test isapprox(csi.a[i] + csi.b[i] + csi.c[i] + csi.d[i], y[i+1], atol=1.e-12)
    end
    # Check meeting knot conditions
    for i = 1:3
        di = csi.b[i+1]
        dip = csi.b[i] + 2*csi.c[i] + 3*csi.d[i]
        @test isapprox(di, dip, atol=1.e-12)
    end
    for i = 1:3
        ddi = 2*csi.c[i+1]
        ddip = 2*csi.c[i]+6*csi.d[i]
        @test isapprox(ddi, ddip, atol=1.e-12)
    end
    # Second derivatives at end points
    @test isapprox(csi.c[1], 0.0, atol = 1.e-12)
    @test isapprox(2*csi.c[4]+6*csi.d[4], 0.0, atol = 1.e-12)
    # Test matching boundary conditions with unequally spaced knots
    x = [0.0, 0.7, 2.3, 3.0, 4.1]
    y = sin.(x)
    csi = CubicSpline(x,y)
    for i = 1:4
        @test csi.a[i] == y[i]
        alpha = x[i+1]-x[i]
        yend = csi.a[i] + csi.b[i]*alpha + csi.c[i]*alpha^2
         + csi.d[i]*alpha^3
        # @test isapprox(yend, y[i+1], atol=1.e-12)
    end
    # Check for continuity near knot 2
    eps = 0.0001
    vl = x[2] - eps
    vg = x[2] + eps
    yl = interp(csi, vl)
    yg = interp(csi, vg)
    @test abs(yl-yg) < 2*eps
    sl = slope(csi, vl)
    sg = slope(csi, vg)
    @test abs(sl-sg) < 2*eps
    sl2 = slope2(csi, vl)
    sg2 = slope2(csi, vg)
    @test abs(sl2-sg2) < 2*eps
    # Check meeting knot conditions
    for i = 1:3
        alpha = x[i+1]-x[i]
        dip = csi.b[i+1]
        di = csi.b[i]+2*csi.c[i]*alpha+3*csi.d[i]*alpha^2
        @test isapprox(di, dip, atol=1.e-12)
    end
    for i = 1:3
        alpha = x[i+1]-x[i]
        ddi = 2*csi.c[i+1]
        ddip = 2*csi.c[i]+6*csi.d[i]*alpha
        @test isapprox(ddi, ddip, atol=1.e-12)
    end
    # Second derivatives at end points
    @test isapprox(csi.c[1], 0.0, atol = 1.e-12)
    alpha = x[5] - x[4]
    @test isapprox(2*csi.c[4]+6*csi.d[4]*alpha, 0.0, atol = 1.e-12)
    end;
end

function graphics_tests()
    x = range(0.0, stop=pi, length=10)
    y = sin.(x)
    cs = CubicSpline(x,y)
    xx = range(0.0, stop=pi, length=97)
    yy = [interp(cs,v) for v in xx]
    yyy = sin.(xx)
    plt.figure()
    plt.plot(x,y,"o",xx,yy,"-",xx,yyy,".")
    plt.title("Regular Interpolation")

    x = cumsum(rand(10));
    x = (x.-x[1]).*pi/(x[10].-x[1])
    y = sin.(x)
    cs = CubicSpline(x,y)
    xx = range(0.0, stop=pi, length=97)
    yy = [interp(cs,v) for v in xx]
    yyy = sin.(xx)
    plt.figure()
    plt.plot(x,y,"o",xx,yy,"-",xx,yyy,".")
    plt.title("Irregular Interpolation, 10 Points")
end

function bump_tests()
    x = [0.0, 0.1, 0.2, 0.3, 0.35, 0.55, 0.65, 0.75];
    y = [0.0, 0.01, 0.02, 0.03, 0.5, 0.51, 0.52, 0.53];
    xx = range(0.0,stop=0.75,length=400);
    sp = CubicSpline(x,y);
    yy = [interp(sp, v) for v in xx]
    pc = pchip(x,y)
    yyy = [interp(pc,v) for v in xx]
    pc2 = pchip2(x,y)
    yyy2 = [interp(pc2,v) for v in xx]
    plt.figure()
    plt.plot(x,y,"o",xx,yy,"-",xx,yyy,"-",xx,yyy2,"-")
    plt.title("Cubic Interpolation")
    plt.legend(("data", "spline", "mean","quad"))
    pc3 = pchip3(x,y)
    yyy3 = [interp(pc3,v) for v in xx]
    plt.figure()
    plt.plot(x, y, "o", xx, yyy3, "-")
    plt.title("PCHIP Interpolation")
    plt.legend(("data", "PCHIP"))
end

function regular_pchip_tests()
    @testset "Regular pchip" begin
    end;
end

function irregular_pchip_tests()
    x=[1.0, 1.8, 2.5, 3.0, 3.9];
    y=cos.(x);
    pc=pchip(x,y)
    @testset "Irregular pchip" begin
    for i=1:5
        # Continuity
        @test interp(pc,x[i])==y[i]
    end
    for i = 2:4
        # Continuity of slope
        eps = 0.000001
        @test isapprox(slope(pc,x[i]-eps),slope(pc,x[i]+eps),atol=4*eps)
    end
    end;
end

if unitTests
    regular_tests()
    irregular_tests()
    regular_pchip_tests()
    irregular_pchip_tests()
end

if graphicsTests
    graphics_tests()
end

if bumpTests
    bump_tests()
end

end # module Interp
