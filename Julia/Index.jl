using Interp
using Test

dotests = true

"""
    Index(n,k)
	
n and k are interpolation objects initialized to index data
The interpolation assumes the wavelength is in Angstroms
"""
struct Index
	n:: Interp.PCHIP
	k:: Interp.PCHIP
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
	Index(pchip(wl,n),pchip(wl,k))
end

"""
    getindex(ndx::Index, wl::Float64)
	
Returns a complex index of refraction at wavelength wl.

* Example
ndx = Index("Al")
n = ndx[431.4]
	
n will have the value of the index of refraction of Al at 431.4 Angstroms 
"""
function Base.getindex(ndx::Index, wl::Float64)
	n = interp(ndx.n, wl)
	k = interp(ndx.k, wl)
	n+k*1im
end

if dotests
	ndx = Index("Al")
	@testset "Index Interpolation" begin
		@test isapprox(interp(ndx.n, 409.8),0.895,atol=0.001)
		@test isapprox(interp(ndx.k, 550.0),2.1e-2,atol=1.e-3)
		nval = ndx[409.8]
		n = real(nval)
		nval = ndx[550.0]
		k = imag(nval)
		@test isapprox(n, 0.895, atol=0.001)
		@test isapprox(k, 2.1e-2, atol=1.e-3)
	end;
end