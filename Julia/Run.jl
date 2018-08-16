using DelimitedFiles

struct Run
	x::Array{Float64,1}
	det::Array{Float64,1}
	m3::Array{Float64,1}
	beam::Array{Float64,1}
	npts::Int
end

function Run(lf ::LogFile, run_number :: Int)
	fname = joinpath(dataDir,filename(lf,run_number))
	data = readdlm(fname, skipstart=1)
	gfact = 10^gain(lf, run_number)
	sr = data[:,2]/gfact
	sm3 = data[:,3]/gfact
	Run(data[:,1],sr,sm3,data[:,4],size(data)[1])
end
