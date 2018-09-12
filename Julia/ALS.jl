module ALS

import Base.length
import PyPlot.plot
export LogFile, Index, Reflectance, Run, filename, comment, runs
export wavelength, rplot, sub, plot2, parratt, gain

dataDir = "X:/ALSData/2016/July2016 ALS data"

include("LogFile.jl")
# Moved Index.jl to Materials in order to standardize and
# localize all of the data io stuff dealing with index
# of refraction in one place.
# include("Index.jl")
include("Run.jl")
include("Reflectance.jl")
include("Refl.jl")
end
