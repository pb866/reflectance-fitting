module ALS

import Base.length
import PyPlot.plot
export LogFile, Index, Reflectance, Run, filename, comment, runs
export wavelength, rplot, sub, plot2, parratt, gain

dataDir = "X:/ALSData/2016/July2016 ALS data"

include("LogFile.jl")
include("Index.jl")
include("Run.jl")
include("Reflectance.jl")
include("Refl.jl")
end
