using DataFrames
using CSV

struct LogFile
	data :: DataFrame
end

# Read in AlF3 data from Zoe's log file (doctored)
function LogFile()
	LogFile(CSV.read(joinpath(dataDir,"July2016.csv"),delim="\t",datarow=33))
end

# Return the comment for a run
function comment(lf ::LogFile, run::Int)
	lf.data[run,:comment]
end

function filename(lf ::LogFile, run::Int)
	lf.data[run,:filename]
end

function gain(lf ::LogFile, run::Int)
	lf.data[run,Symbol("gain/count1")]
end

function x(lf ::LogFile, run::Int)
	lf.data[run, :x]
end

function y(lf ::LogFile, run::Int)
	lf.data[run, :y]
end

function z(lf ::LogFile, run::Int)
	lf.data[run, :z]
end

function runs(lf ::LogFile)
	size(lf.data)[1]
end