"""
Read and interpret CXRO and other data with
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

# types
I've included both abstract and concrete types in this list.

* `OpticalConstant` - general type for optical data
"""
module Materials

using Test
using DataFrames
using CSV
using Interp
using Rakic

export Index
export delta, beta, f1, f2, e2lambda, lambda2e
export delta2, beta2, CXROread, NK, index, Compound

# if true, runs units tests
const dotests = false

abstract type OpticalConstant end;
abstract type DataRange end;
abstract type EnergyOrWavelength <: DataRange end;
abstract type RealIndex <: DataRange end;
abstract type ImagIndex <: DataRange end;
abstract type Nanometer <: EnergyOrWavelength end;
abstract type Angstrom <: EnergyOrWavelength end;
abstract type ElectronVolt <: EnergyOrWavelength end;
abstract type N <: RealIndex end;
abstract type K <: ImagIndex end;
abstract type Delta <: RealIndex end;
abstract type Beta <: ImagIndex end;
abstract type F1 <: RealIndex end;
abstract type F2 <: ImagIndex end;

"""
    Interpolator
# fields
* pchip: (Missing or Interp.PHIP) data to interpolate index data
"""
struct Interpolator
    pchip :: Union{Missing, Interp.PCHIP}
end

struct Index <: OpticalConstant
    ew :: Angstrom
    re :: N
    im :: K
    ninterp :: Interpolator
    kinterp :: Interpolator
end

struct NK <: OpticalConstant
    ew :: Nanometer
    re :: N
    im :: K
    ninterp :: Interpolator
    kinterp :: Interpolator
end

struct F1F2 <: OpticalConstant
    ew :: ElectronVolt
    re :: F1
    im :: F2
    f1interp :: Interpolator
    f2interp :: Interpolator
    A :: Float64
end

struct AngstromRange <: Angstrom
    range::AbstractArray{Float64, 1}
end

struct NanometerRange <: Nanometer
    range::AbstractArray{Float64, 1}
end

struct EnergyRange <: ElectronVolt
    range::AbstractArray{Float64, 1}
end

struct NRange <: N
    range::AbstractArray{Float64,1}
end

struct KRange <: K
    range::AbstractArray{Float64,1}
end

struct F1Range <: F1
    range::AbstractArray{Float64,1}
end

struct F2Range <: F2
    range::AbstractArray{Float64,1}
end

Base.getindex(r :: DataRange, index::Int64) = r.range[index]
Base.length(r::DataRange) = length(r.range)
NK() = NK(EnergyRange([]), NRange([]), KRange([]), Interpolator(missing), Interpolator(missing))
nanometers(r::AngstromRange) = r.range/10.0
nanometers(r::NanometerRange) = r.range
nanometers(r::EnergyRange) = e2lambda.(r.range)
wavelength(oc::OpticalConstant) = nanometers(oc.ew)

"""
    interp(i, x)
interpolate using an `Interpolator` `i` at the point `x`
"""
interp(i::Interpolator, x::Real) = Interp.interp(i.pchip, x)

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
atomic_data["ca"] = AtomicData(1.54,40.078)

const alA=atomic_data["al"].A
const alrho=atomic_data["al"].rho
const fA = atomic_data["f"].A
const frho = atomic_data["f"].rho
const alf3rho = 2.88;

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

"""
    F1F2(element)
Create an F1F2 type by reading the data from a CXRO atomic scattering
    factors data file.

* `element` : lowercase atomic abbreviation of element
"""
function F1F2(element::AbstractString)
    data = CXROread(element)
    ed = Array{Float64}(data[:E])
    f1d = Array{Float64}(data[:f1])
    f2d = Array{Float64}(data[:f2])
    f1interp = pchip3(ed, f1d)
    f2interp = pchip3(ed, f2d)
    F1F2(EnergyRange(data[:E]), F1Range(data[:f1]), F2Range(data[:f2]),
        Interpolator(f1interp), Interpolator(f2interp), atomic_data[element].A)
end

"""
    NK(element, rho, A)
Create an NK type by reading the data from a CXRO atomic scattering
    factors data file.

# Parameters
* `element`: String with the element abbreviation in lower case
* `rho`: Float64 with atomic density in g/cm^3
* `A`: Float64 with atomic mass in g/mole
"""
function NK(element::AbstractString, rho::Float64, A::Float64)
    data = CXROread(element)
    lambda = e2lambda.(data[:E])
    wf = reverse(collect(zip(lambda,data[:f1])))
    n = [1-delta(rho, first(d), A, last(d)) for d in wf]
    ninterp = Interp.pchip3(reverse(lambda),n)
    wf = reverse(collect(zip(lambda,data[:f2])))
    k = [beta(rho, first(d), A, last(d)) for d in wf]
    kinterp = Interp.pchip3(reverse(lambda),k)
    NK(NanometerRange(lambda), NRange(n), KRange(k),
        Interpolator(ninterp), Interpolator(kinterp))
end

function NK(element::AbstractString)
    if haskey(atomic_data,element)
        NK(element, atomic_data[element].rho, atomic_data[element].A)
    else
        error("$element is not in the atomic database")
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

f1(f1f2::F1F2, wl::Float64) = interp(f1f2.f1interp, wl)
f2(f1f2::F1F2, wl::Float64) = interp(f1f2.f2interp, wl)

index(nk::NK, wl::Int64) = getindex(nk, convert(Float64,wl))

function nk_test()
    @testset "NK" begin
    alnk = NK("al")
    @test_throws Exception wnk = NK("w")
    @test isapprox(index(alnk,30.4),0.95123+0.0056009im, atol=1e-5)
    end;
end

function f1f2_test()
    @testset "F1F2" begin
    alf1f2 = F1F2("al")
    @test isapprox(f1(alf1f2,30.4), 2.3066872, atol=1e-7)
    @test isapprox(f2(alf1f2,30.4), 0.1764198, atol=1e-7)
    end;
end

function Compound(rho::Float64, f1f2_1::F1F2, n1::Int64,
        f1f2_2::F1F2, n2::Int64)
    E1 = reverse(f1f2_1.ew.range)
    f1 = reverse(f1f2_1.re.range)
    wl1 = e2lambda.(E1)
    A1 = f1f2_1.A
    wl2 = e2lambda.(reverse(f1f2_2.ew.range))
    wint = pchip3(wl2,reverse(f1f2_2.re.range))
    f2 = [Interp.interp(wint, lambda) for lambda in wl1]
    A2 = f1f2_2.A
    wf = zip(wl1,f1,f2)
    n = [1-delta2(rho, d[1], A1, n1, d[2], A2, n2, d[3]) for d in wf]
    ninterp = pchip3(wl1, n)
    f1 = reverse(f1f2_1.im.range)
    wint = pchip3(wl2,reverse(f1f2_2.im.range))
    f2 = [Interp.interp(wint, lambda) for lambda in wl1]
    wf = zip(wl1,f1,f2)
    k = [beta2(rho, d[1], A1, n1, d[2], A2, n2, d[3]) for d in wf]
    kinterp = pchip3(wl1, k)
    NK(NanometerRange(wl1), NRange(n), KRange(k),
        Interpolator(ninterp), Interpolator(kinterp))
end

function compound_test()
    @testset "compound" begin
    af1f2=F1F2("al")
    ff1f2=F1F2("f")
    cmpd = Compound(2.88, af1f2, 1, ff1f2, 3)
    @test isapprox(index(cmpd,30.4), 0.92994+0.13544im, atol=1e-5)
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

"""
    Index(material::String)

Initialize an Index type from data data in an index file.
"""
function Index(name::String)
	wl, n, k = open(name*".nk") do file
		lines = readlines(file)
		headlines = 1
		# Count header lines to skip them
		while lines[headlines][1] == ';'
			headlines += 1
		end
		headlines -= 1
		# Number of data lines is the rest
		datalines = length(lines) - headlines
		# skip a blank line at the end if there is one
		if length(lines[end]) == 0
			datalines -= 1
		end
		wl = zeros(datalines)
		n = zeros(datalines)
		k = zeros(datalines)
		for i in 1:datalines
			sline = lines[i+headlines]
			tokens = split(sline)
			wlstr = tokens[1]
			nstr = tokens[2]
			kstr = tokens[3]
			wl[i] = parse(Float64, wlstr)
			n[i] = parse(Float64, nstr)
			k[i] = parse(Float64, kstr)
		end
		(wl, n, k)
	end
	Index(AngstromRange(wl), NRange(n), KRange(k), Interpolator(pchip3(wl,n)),
        Interpolator(pchip3(wl,k)))
end

"""
    getindex(ndx::Index, wl::Float64)

Returns a complex index of refraction at wavelength wl nm.

* Example
ndx = Index("Al")
n = ndx[431.4]

n will have the value of the index of refraction of Al at 431.4 Angstroms
"""
function Base.getindex(ndx::Index, wl::Float64)
    # data has index in Angstroms, convert wl to Angstroms before interp
	n = interp(ndx.ninterp, wl*10.0)
	k = interp(ndx.kinterp, wl*10.0)
	n+k*1im
end

Base.getindex(ndx::Index, wl::Int64) = getindex(ndx, convert(Float64,wl))

function index_tests()
	ndx = Index("Al")
	@testset "Index Interpolation" begin
		@test isapprox(interp(ndx.ninterp, 409.8),0.895,atol=0.001)
		@test isapprox(interp(ndx.kinterp, 550.0),2.1e-2,atol=1.e-3)
		nval = ndx[40.98]
		n = real(nval)
		nval = ndx[55.00]
		k = imag(nval)
		@test isapprox(n, 0.895, atol=0.001)
		@test isapprox(k, 2.1e-2, atol=1.e-3)
	end;
end

function wavelength_test()
    # Create some fake data to test
    imiss = Interpolator(missing)
    ndx = Index(AngstromRange([250,260]), NRange([0.8,0.9]), KRange([0.1,0.2]),
        imiss, imiss)
    nk = NK(NanometerRange([18,19]), NRange([0.8,0.9]), KRange([0.1,0.2]),
        imiss, imiss)
    f1f2 = F1F2(EnergyRange([21.5,32.7]), F1Range([0.8, 1.2]),
        F2Range([2.5, 0.5]), imiss, imiss, 2.5)
    @testset "wavelength" begin
    @test wavelength(ndx)[1] == 25
    @test wavelength(nk)[2] == 19
    @test isapprox(wavelength(f1f2)[1], 1239.8/21.5, rtol = 5e-5)
    end;
end

if dotests
    @testset "Materials" begin
    wavelength_test()
    index_tests()
    conversion_test()
    molecular_conversion_test()
    read_test()
    nk_test()
    f1f2_test()
    compound_test()
    end;
end

end
