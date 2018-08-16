module ALS

export LogFile, Index, Reflectance, Run, filename, gain, x, y, z, comment, runs, wavelength

dataDir = "X:/ALSData/2016/July2016 ALS data"

include("LogFile.jl")
include("Index.jl")
include("Reflectance.jl")
include("Run.jl")

end