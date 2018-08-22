module Interp
# Code for interpolation for various orders
using LinearAlgebra
using Test
using PyPlot
import Base.length

unitTests = true
graphicsTests = false

"""
    CubicSpline(x,a,b,c,d)

concrete type for holding the data needed
    to do a cubic spline interpolation
"""
struct CubicSpline
    x
    a::Array{Float64,1}
    b::Array{Float64,1}
    c::Array{Float64,1}
    d::Array{Float64,1}
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

"""
    interp(cs::CubicSpline, v::Float)

Interpolate to the value corresonding to v

# Examples
```
x = cumsum(rand(10))
y = cos.(x);
cs = CubicSpline(x,y)
v = interp(ds, 1.2)
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

function do_tests()
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
        @test isapprox(yend, y[i+1], atol=1.e-12)
    end
    # Check meeting knot conditions
    # All of the test in the following loop fail
    for i = 1:4
        alpha = x[i+1]-x[i]
        dip = csi.b[i+1]
        di = csi.b[i]+2*cs.c[i]*alpha+3*cs.d[i]*alpha^2
        @test isapprox(di, dip, atol=1.e-12)
        ddi = 2*csi.c[i+1]
        ddip = 2*csi.c[i]+6*csi.d[i]*alpha
        @test isapprox(ddi, ddip, atol=1.e-12)
    end
    # Second derivatives at end points
    @test isapprox(csi.c[1], 0.0, atol = 1.e-12)
    @test isapprox(2*csi.c[4]+6*csi.d[4], 0.0, atol = 1.e-12)
    end;
end

function graphics_tests()
    x = range(0.0, stop=pi, length=10)
    y = sin.(x)
    cs = CubicSpline(x,y)
    xx = range(0.0, stop=pi, length=97)
    yy = [interp(cs,v) for v in xx]
    yyy = sin.(xx)
    figure()
    plot(x,y,"o",xx,yy,"-",xx,yyy,".")
    title("Regular Interpolation")

    x = cumsum(rand(10));
    x = (x.-x[1]).*pi/(x[10].-x[1])
    y = sin.(x)
    cs = CubicSpline(x,y)
    xx = range(0.0, stop=pi, length=97)
    yy = [interp(cs,v) for v in xx]
    yyy = sin.(xx)
    figure()
    plot(x,y,"o",xx,yy,"-",xx,yyy,".")
    title("Irregular Interpolation, 10 Points")

    x = cumsum(rand(40));
    x = (x.-x[1]).*pi/(x[40].-x[1])
    y = sin.(x)
    cs = CubicSpline(x,y)
    xx = range(0.0, stop=pi, length=397)
    yy = [interp(cs,v) for v in xx]
    yyy = sin.(xx)
    figure()
    plot(x,y,"o",xx,yy,"-",xx,yyy,".")
    title("Irregular Interpolation, 40 Points")
end

struct Hermite
    x::Array{Float64, 1}
    y::Array{Float64,1}
    d::Array{Float64,1}
end

function Hermite(x::Array{Float64,1}, y::Array{Float64,1})
    d=zeros(length(x))
    # d[1] =
end

if unitTests
    do_tests()
end

if graphicsTests
    graphics_tests()
end

end # module Interp
