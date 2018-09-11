# Rakic calculations
module Rakic

using Interp
using Test

const dotests = false

export alindex

function CSVread(filename)
	f = open(filename)
	lines = readlines(f)
	npts = length(lines)-1
	eV = Array{Float64}(undef, npts)
	n = Array{Float64}(undef, npts)
	k = Array{Float64}(undef, npts)
	for i=1:npts
		parts = split(lines[i+1],",")
		eV[i] = parse(Float64, parts[1])
		n[i] = parse(Float64, parts[2])
		k[i] = parse(Float64, parts[3])
	end
	(eV, n, k)
end

data_path = homedir()*"/Documents/alf3-paper/papers/Rakic.csv"
(eV, n, k) = CSVread(data_path)
ninterp = Interp.pchip3(eV, n)
kinterp = Interp.pchip3(eV, k)

hc = 1239.84197
lambda2e(lambda) = hc/lambda

function alindex(wl::Float64)
    eV = lambda2e(wl)
    n = interp(ninterp, eV)
    k = interp(kinterp, eV)
    n+k*1im
end

if dotests
    @testset "Rakic" begin
    nk = alindex(30.4)
    @test isapprox(nk,0.94492866+0.007637106im,atol=1e-6)
    end;
end

end
