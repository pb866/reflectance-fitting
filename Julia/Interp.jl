module Interp
# Code for interpolation for various orders
using LinearAlgebra
using Test

unit_tests = true

"""
    CubicSpline(x,a,b,c,d)

concrete type for holding the data needed
    to do a cubic spline interpolation
"""
struct CubicSpline
    x::Array{Float64,1}
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
    yy[1] = 3*y[2]-y[1]
    du[1] = 1
    dd[1] = 2
    for i=2:len-1
        yy[i] = 3*(y[i+1]/alpha[i]^2+
            y[i]*(alpha[i-1]^(-2)+alpha[i]^(-2))-
            y[i-1]/alpha[i-1]^2)
        dl[i-1] = alpha[i-1]^(-2)
        du[i] = dl[i-1]
        dd[i] = 2*(alpha[i-1]^(-2)+alpha[i]^(-2))
    end
    yy[len] = 3*(y[len]-y[len-1])
    dl[len-1] = 1
    dd[len] = 2
    # Solve the tridiagonal system for the derivatives D
    dm = Tridiagonal(dl,dd,du)
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
    t = (v-cs.x[segment])/(cs.x[segment+1]-cs.x[segment])
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

function do_tests()
    @testset "interpolation" begin
    # Test not enough points exception
    x = [1.0, 2.0]
    y = [2.0, 4.0]
    @test_throws ErrorException CubicSpline(x,y)
    x = [0.2, 1.4, 3.8, 5.7]
    y = [1.5, 3.0, 3.7, 2.5]
    cs = CubicSpline(x,y)
    @test_throws ErrorException interp(cs, 0.0)
    @test_throws ErrorException interp(cs, 6.0)
    @test region(x, 0.3) == 1
    @test region(x, 0.2) == 1
    @test region(x, 5.7) == 3
    @test region(x, 2.1) == 2
    @test region(x, 4.0) == 3
    @test interp(cs, 0.2) == 1.5
    @test interp(cs, 1.4) == 3.0
    @test isapprox(interp(cs, 5.7), 2.5, atol=1e-14)
    end;
end

if unit_tests
    do_tests()
end

end # module Interp
