# for Debugging Reflectance constructor
# and plot function

using ALS
using Serialization

# lf = LogFile()
f = open("logfile.bin")
lf = deserialize(f)
close(f)

irun = Run(lf, 20)
i0run = Run(lf, 11)
darkrun = Run(lf,12)
dark0run = Run(lf,13)
lambda = wavelength(lf, 20)

rfl = Reflectance(irun, i0run, darkrun, dark0run, lambda)

# ALS.plot(rfl)

# For sample A I don't need this function. I will wait to fully
# develop it until later.

function +(r1::Reflectance, r2::Reflectance)
    if(abs(r1.wavelength-r2.wavelength)>0.1)
        error("Wavelengths don't match.")
    end
    Reflectance(vcat(r1.theta, r2.theta), vcat(r1.refl, r2.refl), r1.wavelength)
end
