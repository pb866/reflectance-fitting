"""
Read and interpret CXRO data with
various units.

# functions:
* `delta` - compute delta from f1
* `beta` - compute beta from f2
* `f1` - compute f1 from delta
* `f2` - compute f2 from beta
* `e2lambda` - convert from energy to wavelength
* `lambda2e` - convert from wavelength to energy
* `delta2` - delta for a binary compound
* `beta2` - beta from a binary compound
"""
module CXRO

using Test
using DataFrames
using CSV
using Interp
using Rakic

# if true, runs units tests
const dotests = false
const alA=26.982
const alrho=2.70
const fA = 18.998
const frho = 1.0 # dummy value, it is a gas...
const alf3rho = 2.88;

export delta, beta, f1, f2, e2lambda, lambda2e
export delta2, beta2, CXROread, NK, index, Compound

"""
    delta(rho,lambda,A,f1)
compute delta from f1
# Parameters
* `rho` : density in g/cm^2
* `lambda` : wavelength in nm
* `A` : atomic mass in g/mole
* `f1` : real part of oscillator strength
"""
delta(rho, lambda, A, f1) = 2.7008644e-4*rho*lambda^2/A*f1

"""
    beta(rho,lambda,A,f2)
compute beta from f2
# Parameters
* `rho` : density in g/cm^2
* `lambda` : wavelength in nm
* `A` : atomic mass in g/mole
* `f2` : imaginary part of oscillator strength
"""
beta(rho, lambda, A, f2) = 2.7008644e-4*rho*lambda^2/A*f2

"""
    f1(rho,lambda,A,delta)
compute f1 from delta
# Parameters
* `rho` : density in g/cm^2
* `lambda` : wavelength in nm
* `A` : atomic mass in g/mole
* `delta` : real part of index of refraction
"""
f1(rho, lambda, A, delta) = 3702.518*A/(rho*lambda^2)*delta

"""
    f2(rho,lambda,A,beta)
compute f1 from delta
# Parameters
* `rho` : density in g/cm^2
* `lambda` : wavelength in nm
* `A` : atomic mass in g/mole
* `beta` : imaginary part of index of refraction
"""
f2(rho, lambda, A, beta) = 3702.518*A/(rho*lambda^2)*beta

hc = 1239.84197
"""
    e2lambda(e)
convert from eV to nm

`e` is the energy in eV

function returns wavelength in nm
"""
e2lambda(e) = hc/e

"""
    lambda2e(lambda)
convert from eV to nm

`lambda` is the wavelength in nm

function returns enery in eV
"""
lambda2e(lambda) = hc/lambda

function conversion_test()
    @testset "conversion" begin
    # Tests for Al at 30.2555 eV
    alf1=2.31139 # f2 for Al at 30.2555 eV
    alf2=0.175442 # f2 for Al at 30.2555 eV
    alA=26.982 # atomic mass of Al
    alrho=2.70 # density of Al
    E=30.2555 # energy for this test
    lambda=e2lambda(E)
    @test isapprox(lambda, 40.97906, atol = 1e-5)
    d = delta(alrho, lambda, alA, alf1)
    b = beta(alrho, lambda, alA, alf2)
    @test isapprox(d,0.1049094, atol = 1e-5)
    @test isapprox(b, 0.007962753, atol = 1e-6)
    @test isapprox(f1(alrho, lambda, alA, d), alf1, atol = 1e-6)
    @test isapprox(f2(alrho, lambda, alA, b), alf2, atol = 1e-6)
    end;
end

"""
    delta2(rho, lambda, A1, n1, f11, A2, n2, f12)

compute the real part of the index of refaction for a
    binary compound consisting of `n1` atoms of
    the first element and `n2` atoms of the second
    element.
# Parameters
* `rho` - density of the compound in g/cm^3
* `lambda` - wavelength in nm
* `A1` - atomic mass of first element in g/mole
* `n1` - number of atoms of first element
* `f11` - f1 for first element
* `A2` - atomic mass of second element in g/mole
* `n2` - number of atoms of second element in g/mole
* `f12` - f1 for second element
# Example
Computing delta for AlF3
```
alf1=2.31139 # f1 for Al
alA=26.982   # A for Al
ff1=0.517012 # f1 for F
fA=18.998    # A for F
alf3rho=2.88 # density of AlF3
E=30.2555    # energy in eV
lambda=e2lambda(E) # wavelength in nm
delta2(alf3rho, lambda, alA, 1, alf1, fA, 3, ff1)
```
returns 0.060079215657434994
"""
function delta2(rho, lambda, A1, n1, f11, A2, n2, f12)
    amol = n1*A1+n2*A2
    d1 = delta(rho, lambda, amol, f11) * n1
    d2 = delta(rho, lambda, amol, f12) * n2
    d1+d2
end

"""
    beta2(rho, lambda, A1, n1, f21, A2, n2, f22)

compute the the imaginary part of
    the index of refaction for a
    binary compound consisting of `n1` atoms of
    the first element and `n2` atoms of the second
    element.
# Parameters
* `rho` - density of the compound in g/cm^3
* `lambda` - wavelength in nm
* `A1` - atomic mass of first element in g/mole
* `n1` - number of atoms of first element
* `f21` - f2 for first element
* `A2` - atomic mass of second element in g/mole
* `n2` - number of atoms of second element in g/mole
* `f22` - f2 for second element
# Example
Computing beta for AlF3
```
alf2=0.175442 # f2 for Al
alA=26.982    # A for Al
ff2=4.96088   # f2 for F
fA=18.998     # A for F
alf3rho=2.88  # density of AlF3
E=30.2555     # energy in eV
lambda=e2lambda(E) # wavelength in nm
beta2(alf3rho, lambda, alA, 1, alf2, fA, 3, ff2)
```
returns 0.23422526564012883
"""
function beta2(rho, lambda, A1, n1, f21, A2, n2, f22)
    amol = n1*A1+n2*A2
    b1 = beta(rho, lambda, amol, f21) * n1
    b2 = beta(rho, lambda, amol, f22) * n2
    b1+b2
end

function molecular_conversion_test()
    @testset "molecular conversion" begin
    alf1=2.31139
    alf2=0.175442
    alA=26.982
    ff1=0.517012
    ff2=4.96088
    fA=18.998
    alf3rho=2.88
    E=30.2555
    lambda=e2lambda(E)
    d = delta2(alf3rho, lambda, alA, 1, alf1, fA, 3, ff1)
    b = beta2(alf3rho, lambda, alA, 1, alf2, fA, 3, ff2)
    @test isapprox(d,0.06007458, atol = 5e-6)
    @test isapprox(b,0.2342305, atol = 6e-6)
    end;
end

"""
    CXROread(element)
Read f1 and f2 data from CXRO data file stored
in homedir()/Documents/alf3-paper/analysis/CXRO.
If the data file isn't stored locally, it is downloaded
from the Henke database at CXRO.

`element` is the lower case string with the abbreviation
of the element.

returns a DataFrame with columns `E`, `f1`, and `f2`.

`E` is the energy in eV
# Example
```
aldf = CXROread("al");
aldf[1,:]
1×3 DataFrames.DataFrame
│ Row │ E    │ f1      │ f2     │
├─────┼──────┼─────────┼────────┤
│ 1   │ 10.0 │ -9999.0 │ 3.1199 │
```
"""
function CXROread(element::AbstractString)
    data_dir = homedir()*"/Documents/alf3-paper/analysis/CXRO/"
    filename = data_dir*element*".nff"
    if !isfile(filename)
        henke_url = "http://henke.lbl.gov/optical_constants/sf/"
        urlname = henke_url*element*".nff"
        download(urlname, filename)
    end
    lines = open(data_dir*element*".nff") do f
        readlines(f)
    end
    npts = length(lines)-1
    E = Array{Union{Missing, Float64}}(missing,npts)
    f1 = Array{Union{Missing, Float64}}(missing,npts)
    f2 = Array{Union{Missing, Float64}}(missing,npts)
    for i = 1:npts
        parts = split(lines[i+1])
        E[i] = parse(Float64,parts[1])
        f1[i] = parse(Float64,parts[2])
        f2[i] = parse(Float64, parts[3])
    end
    DataFrame(E=E,f1=f1,f2=f2)
end

function read_test()
    @testset "read CXRO data" begin
    al_df = CXROread("al")
    @test al_df[1,:E] == 10.0
    @test al_df[1,:f1] == -9999.0
    @test al_df[1,:f2] == 3.1199
    @test length(al_df) == 3
    @test size(al_df,1) == 504
    end;
end

# Atomic database
struct AtomicData
    rho :: Float64
    A :: Float64
end
atomic_data = Dict{AbstractString,AtomicData}()
atomic_data["al"] = AtomicData(2.70,26.982)
atomic_data["si"] = AtomicData(2.328,28.09)
atomic_data["f"] = AtomicData(-1,18.998)
atomic_data["o"] = AtomicData(-1,15.9994)

struct NK
    data :: DataFrame
    ninterp :: Union{Missing, Interp.PCHIP}
    kinterp :: Union{Missing, Interp.PCHIP}
    A :: Union{Missing, Float64}
end

function NK(element::AbstractString, rho::Float64, A::Float64)
    data = CXROread(element)
    lambda = e2lambda.(data[:E])
    wf = reverse(collect(zip(lambda,data[:f1])))
    n = [1-delta(rho, first(d), A, last(d)) for d in wf]
    ninterp = Interp.pchip3(reverse(lambda),n)
    wf = reverse(collect(zip(lambda,data[:f2])))
    k = [beta(rho, first(d), A, last(d)) for d in wf]
    kinterp = Interp.pchip3(reverse(lambda),k)
    NK(data, ninterp, kinterp, A)
end

function NK(element::AbstractString)
    if haskey(atomic_data,element)
        if atomic_data[element].rho > 0
            NK(element, atomic_data[element].rho,
                atomic_data[element].A)
        else
            NK(CXROread(element), missing, missing,
                atomic_data[element].A )
        end
    else
        NK(CXROread(element), missing, missing, missing)
    end
end

"""
    index(nk::NK, wl::Float64)

Returns a complex index of refraction at wavelength wl nm.

* Example
nk = NK("al")
n = index(nk, 30.4)

n will have the value of the index of refraction of Al at 30.4 nm
(0.95123+0.056009im)
"""
function index(nk::NK, wl::Float64)
    # data has index in eV, convert wl to eV before interp
	n = interp(nk.ninterp, wl)
	k = interp(nk.kinterp, wl)
	n+k*1im
end

function nk_test()
    @testset "NK" begin
    alnk = NK("al")
    wnk = NK("w")
    @test alnk.A == 26.982
    @test ismissing(wnk.ninterp)
    @test ismissing(wnk.kinterp)
    @test ismissing(wnk.A)
    @test isapprox(index(alnk,30.4),0.95123+0.0056009im, atol=1e-5)
    end;
end

struct Compound
    ninterp :: Interp.PCHIP
    kinterp :: Interp.PCHIP
end

function Compound(rho::Float64, nk1::NK, n1::Int64,
        nk2::NK, n2::Int64)
    E1 = reverse(nk1.data[:E])
    f1 = reverse(nk1.data[:f1])
    wl1 = e2lambda.(E1)
    A1 = nk1.A
    pchip(x::Array{Float64},
        y::Array{Union{Missing,Float64}}) =
        pchip3(Array{Float64}(x), Array{Float64}(y))
    wl2 = e2lambda.(reverse(nk2.data[:E]))
    wint = pchip(wl2,reverse(nk2.data[:f1]))
    f2 = [interp(wint, lambda) for lambda in wl1]
    A2 = nk2.A
    wf = zip(wl1,f1,f2)
    n = [1-delta2(rho, d[1], A1, n1, d[2], A2, n2, d[3]) for d in wf]
    ninterp = pchip3(wl1, n)
    f1 = reverse(nk1.data[:f2])
    wint = pchip(wl2,reverse(nk2.data[:f2]))
    f2 = [interp(wint, lambda) for lambda in wl1]
    wf = zip(wl1,f1,f2)
    k = [beta2(rho, d[1], A1, n1, d[2], A2, n2, d[3]) for d in wf]
    kinterp = pchip3(wl1, k)
    Compound(ninterp, kinterp)
end

function index(nk::Compound, wl::Float64)
    # data has index in eV, convert wl to eV before interp
	n = interp(nk.ninterp, wl)
	k = interp(nk.kinterp, wl)
	n+k*1im
end

index(nk::Union{NK, Compound}, wl::Int64) = getindex(nk, convert(Float64,wl))

function compound_test()
    @testset "compound" begin
    ank=NK("al")
    fnk=NK("f")
    cmpd = Compound(2.88, ank, 1, fnk, 3)
    @test isapprox(index(cmpd,30.4), 0.92994+0.13544im, atol=1e-5)
    end;
end

if dotests
    @testset "CXRO" begin
    conversion_test()
    molecular_conversion_test()
    read_test()
    nk_test()
    compound_test()
    end;
end

function rakic()
    A = 26.986
    rho = 2.70
    af1 = [f1(rho, e2lambda(E), A, 1 .- Rakic.n) for E in Rakic.eV]
    af2 = [f2(rho, lambda, A, Rakic.k) for E in Rakic.eV]
    df = DataFrame(E=Rakic.eV, af1, af2)
    lambda = e2lambda.(reverse(Rakic.eV))
    rn = reverse(Rakic.n)
    rk = reverse(Rakic.k)
    ninterp = Interp.pchip3(lambda, rn)
    kinterp = Interp.pchip3(lambda, rk)
    NK(df, rn, rk, A)
end

end
