using Interp
using Test

# If true, run unit tests
unit_tests = false

"""
    parratt(n, x, thetad, lam, fractions, sigma)

Reflectance from a multilayer mirror

# PARAMETERS
* n : array of index of refractions for stack starting with index of the incident layer
* x : array of layer thicknesses in nm for stack starting with the thickness of the incident layer (usually vacuum). The vacuum and substrate thicknesses should be set to 0.
* thetad::Float64
    incident angle in degrees
* lam ::Float64
    wavelength in nanometers
* fractions::Float64
    fraction of light with s polarization. If zero, this is calculated
    using Gullikson formula for synchrotron
* sigma:: array of number
    interface roughness in nm

# RETURNS
* ref: number
        reflectance of the layer for s polarization

# EXAMPLE
```
lam = 15.0
AlIndex = Index("Al")
alndx = AlIndex[lam]
SiO2Index = Index("SiO2")
sio2ndx = SiO2Index[lam]
n = [1.0, alndx, sio2ndx]
t = [0.0, 20.0, 0.0]
thetad = 20.0
p = parratt(n, t, thetad, lam)
0.012606490280013215
```
"""
function parratt(n::Array{Complex{Float64},1},
    x::Array{Float64,1},
    thetad::Float64, lam::Float64; fractions::Float64=0.0,
    sigma::Union{Array{Float64,1},Float64}=0.0)
    if fractions==0
        fractions = fracs(lam)
    end
    fractionp = 1-fractions
    S = sqrt.(n.^2 .- cos(thetad*pi/180)^2)
    k = 2*pi/lam
    C = exp.(2im.*S.*x.*k)
    rs = 0.0
    rp = 0.0

    qz = k*sin(thetad*pi/180) # for Debye-Waller correction
    if typeof(sigma) == Float64
        sigma = ones(length(n)-1)*sigma
    end
    eta = exp.(-2.0*qz^2 .* sigma.^2) # Debye-Waller roughness correction

    for m in range(length(n), stop=2, step=-1)
        fs = (S[m-1]-S[m])/(S[m-1]+S[m])
        fp = ((n[m]^2 * S[m-1] - n[m-1]^2 * S[m])/
              (n[m]^2 * S[m-1] + n[m-1]^2 * S[m]))
        rs = C[m-1]*(fs*eta[m-1]+rs*eta[m-1]^2)/(1.0+fs*rs*eta[m-1])
        rp = C[m-1]*(fp*eta[m-1]+rp*eta[m-1]^2)/(1.0+fp*rp*eta[m-1])
    end
    fractionp*abs2(rp)+fractions*abs2(rs)

end

function test_parratt()
    @testset "Parratt" begin
    lam = 15.0
    AlIndex = Index("Al")
    alndx = AlIndex[lam]
    SiO2Index = Index("SiO2")
    sio2ndx = SiO2Index[lam]
    n = [1.0, alndx, sio2ndx]
    t = [0.0, 20.0, 0.0]
    thetad = 20.0
    p = parratt(n, t, thetad, lam)
    @test isapprox(p, 0.012606490280013215, atol=1e-5)
    end;
end

ipts = [1.05685E+1 9.37888E-1;
     1.40446E+1 9.26536E-1;
     1.75270E+1 9.15094E-1;
     2.25761E+1 9.02016E-1;
     2.95528E+1 8.85600E-1;
     3.80743E+1 8.70842E-1;
     4.67841E+1 8.57696E-1;
     5.74863E+1 8.44550E-1;
     7.06822E+1 8.26362E-1;
     8.69072E+1 8.08173E-1;
     1.15690E+2 7.83378E-1;
     1.46820E+2 7.63555E-1;
     1.83342E+2 7.47070E-1;
     2.40104E+2 7.27294E-1;
     3.04907E+2 7.02429E-1;
     3.99133E+2 6.86013E-1;
     5.14002E+2 6.74616E-1;
     6.72413E+2 6.63242E-1;
     9.06564E+2 6.61998E-1;
     1.04406E+3 6.63885E-1];

 xfr = ipts[:,1]
 yfr = ipts[:,2]
 pcs = pchip(log10.(xfr),yfr)
 """
    fracs(lam::Float64)
 fraction of s poplarization at the ALS

# Parameters
 lam : number
     wavelength in nm

# Returns
 fr : number
     fraction of s polarization at this wavelength

# Example
```
>>>fracs(25.1)
0.9272135287221943
```
"""
function fracs(lam::Float64)
     eV=log10(1239.8/lam)
     return (interp(pcs,eV)+1)/2
end

function test_fracs()
    @testset "fracs" begin
    @test isapprox(fracs(25.1), 0.92717938, atol=1e-7)
    end;
end

if unit_tests
    test_fracs()
    test_parratt()
end
