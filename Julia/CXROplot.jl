# Plot Al and AlF3 data using CXRO database

using Materials
import PyPlot
using LaTeXStrings
import ALS
import CSV
import Interp

const plt = PyPlot

function aplot()
    ank = NK("al")
    wl = 6.3:0.1:30.4
    en = [index(ank, lambda) for lambda in wl]
    n = real.(en)
    k = imag.(en)
    plt.figure()
    plt.plot(wl,n)
    plt.title("Aluminum Index of Refraction from CXRO")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("n")
    plt.figure()
    plt.plot(wl,k)
    plt.title("Aluminum Index of Refraction from CXRO")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("k")
end

function alf3plot()
    ank = NK("al")
    fnk = NK("f")
    rho = 2.88 # density of AlF3
    alf3 = Compound(rho, ank, 1, fnk, 3)
    wl = 6.3:0.1:30.4
    en = [index(alf3, lambda) for lambda in wl]
    n = real.(en)
    k = imag.(en)
    plt.figure()
    plt.plot(wl,n)
    plt.title("Alumina Index of Refraction from CXRO")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("n")
    plt.figure()
    plt.plot(wl,k)
    plt.title("Alumina Index of Refraction from CXRO")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("k")
end

function volta()
    andx = ALS.Index("Al")
    ank = NK("al")
    wl = 6.3:0.1:80
    en = [index(ank, lambda) for lambda in wl]
    ven = [andx[lambda] for lambda in wl]
    n = real.(en)
    k = imag.(en)
    vn = real.(ven)
    vk = imag.(ven)
    # Restrict real part to where there is valid CXRO data
    maxwl = 6.3
    wln = zip(wl,n)
    for (lambda, n) in wln
        if n<10.0 maxwl = lambda end
    end
    println("Maximum wavelength for real index is $maxwl")
    wl2 =[lambda for lambda in 6.3:0.1:maxwl]
    en2 = [index(ank, lambda) for lambda in wl2]
    ven2 = [andx[lambda] for lambda in wl2]
    n2 = real.(en2)
    vn2 = real.(ven2)
    plt.figure()
    plt.plot(wl2,n2,wl2,vn2)
    plt.title("Aluminum Index of Refraction")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("n")
    plt.legend(["CXRO","volta"])
    plt.figure()
    plt.plot(wl,k,wl,vk)
    plt.title("Aluminum Index of Refraction")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("k")
    plt.legend(["CXRO","volta"])
end

function rakic()
	data_path = homedir()*"/Documents/alf3-paper/papers/Rakic.csv"
	rdf = CSV.read(data_path)
	pchip(x::Array{Union{Missing,Float64}}, y::Array{Union{Missing,Float64}}) =
		Interp.pchip3(Array{Float64}(x), Array{Float64}(y))
	ninterp = pchip(rdf[:eV], rdf[:n])
	kinterp = pchip(rdf[:eV], rdf[:k])

	hc = 1239.84197
	lambda2e(lambda) = hc/lambda

	function alindex(wl::Float64)
		eV = lambda2e(wl)
		n = Interp.interp(ninterp, eV)
		k = Interp.interp(kinterp, eV)
		n+k*1im
	end

    ank = NK("al")
    wl = 6.3:0.1:80
    en = [index(ank, lambda) for lambda in wl]
    n = real.(en)
    k = imag.(en)
    # Restrict real part to where there is valid CXRO data
    maxwl = 6.3
    wln = zip(wl,n)
    for (lambda, n) in wln
        if n<10.0 maxwl = lambda end
    end
    println("Maximum wavelength for real index is $maxwl")
    wl2 =[lambda for lambda in 6.3:0.1:maxwl]
    en2 = [index(ank, lambda) for lambda in wl2]
    n2 = real.(en2)
    rn = [alindex(lambda) for lambda in wl]
    rn2 = real.(rn)
    plt.figure()
    plt.plot(wl2,n2,wl,rn)
    plt.title("Aluminum Index of Refraction")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("n")
    plt.legend(["CXRO","Rakic"])
    en = [index(ank, lambda) for lambda in wl]
    ck = imag.(en)
    rn = [alindex(lambda) for lambda in wl]
    rk = imag.(rn)
    plt.figure()
    plt.plot(wl,ck,wl,rk)
    plt.title("Aluminum Index of Refraction")
    plt.xlabel(L"$\lambda$, nm")
    plt.ylabel("k")
    plt.legend(["CXRO","Rakic"])
end

# aplot()
# alf3plot()
# volta()
rakic()
