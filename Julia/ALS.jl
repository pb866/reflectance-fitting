module ALS

import Base.length
import PyPlot.plot
export LogFile, Index, Reflectance, Run, filename, gain, x, y, z, comment, runs
export wavelength, plot, sub, plot2

dataDir = "X:/ALSData/2016/July2016 ALS data"

include("LogFile.jl")
include("Index.jl")
include("Run.jl")
include("Reflectance.jl")
include("Refl.jl")
end
