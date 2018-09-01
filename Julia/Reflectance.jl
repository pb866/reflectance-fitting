using Interp
import PyPlot
import Base.length
using LaTeXStrings

"""
Compound type with three fields:

  * theta: angles at which reflectance is known
  * refl : reflectance at those angles
  * wavelength: wavelength for this reflectnace

The following methods are defined to work with of this type:
    length, sub, +, Reflectance, plot
"""
struct Reflectance
  theta::Array{Float64,1}
  refl::Array{Float64,1}
  wavelength::Float64
end

length(r::Reflectance) = length(r.theta)

"""
Construct a Reflectance type from Run data
"""
function Reflectance(irun::Run, i0run::Run, darkrun::Run,
    dark0run::Run, lambda::Float64)
  idata = irun.det
  theta = irun.x
  i0wavelength = i0run.x
  i0data = i0run.det
  i0_interp = pchip(i0wavelength, i0data)
  i0 = interp(i0_interp,lambda)
  mean(run::Run) = sum(run.det)/length(run.det)
  dark = mean(darkrun)
  dark0 = mean(dark0run)
  refl = (idata.-dark)./(i0-dark0)
  Reflectance(theta, refl, lambda)
end

function rplot(refl::Reflectance)
    PyPlot.figure()
    PyPlot.semilogy(refl.theta, refl.refl,"o")
    wl = round(refl.wavelength, digits=2)
    PyPlot.title("AlF\$_3\$ Reflectance at $wl nm")
    PyPlot.xlabel(L"$\theta$, degrees")
    PyPlot.ylabel("reflectance")
end

function plot2(refl::Reflectance, reflb::Reflectance)
    PyPlot.figure()
    PyPlot.semilogy(refl.theta, refl.refl,"ro", reflb.theta, reflb.refl, "go")
    wl = round(refl.wavelength, digits=2)
    PyPlot.title("AlF\$_3\$ Reflectance at $wl nm")
    PyPlot.xlabel("angle, degrees")
    PyPlot.ylabel("reflectance")
end

"""
Create a subset of the points in Reflectance r.

Syntax: sub(r::Reflectance; min::Float64=-1.0, max::Float64=400.0)

The subset has data for theta between min and max inclusive.
"""
function sub(r::Reflectance; min::Float64=-1.0, max::Float64=400.0)
    th = Float64[]
    rf = Float64[]
    for i = 1:length(r.theta)
        if (r.theta[i]>=min) & (r.theta[i] <= max)
            th = cat(th, r.theta[i], dims=1)
            rf = cat(rf, r.refl[i], dims=1)
        end
    end
    Reflectance(th, rf, r.wavelength)
end
