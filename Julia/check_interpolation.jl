# process AlF3 data
# need to go to this directory first with als()
using ALS
using Serialization
using PyPlot
using Interpolations

# If true, the data will be reinitialized
initRead = false

if initRead
  # To initially read everything
  # This will take a long time
  lf = LogFile()
  alldata = [Run(lf, rnum) for rnum = 1:ALS.runs(lf)];
  # Serialize it so it is faster next time
  f=open("alldata.bin",write=true)
  serialize(f, alldata)
  close(f)
  f=open("logfile.bin",write=true)
  serialize(f, lf)
  close(f)
else
  # Faster way if you already have the data
  f = open("alldata.bin")
  alldata = deserialize(f)
  close(f)
  f = open("logfile.bin")
  lf = deserialize(f)
  close(f)
end

# The first I0 run is 11
w11=alldata[11].x
r11=alldata[11].det
# Compare w11 with a range
nel = size(w11,1)
rw11=collect(range(w11[1], stop=w11[nel], length=nel));
dif11=w11.-rw11
figure()
plot(w11, dif11,"o")
title("Difference between Data and Range for Run 11")
xlabel("wavelength, nm")
ylabel("difference, nm")
figure()
dw = w11[50]-49.5
rw11=collect(range(25.0+dw, stop=49.5+dw, length=nel))
dif11=w11.-rw11
plot(w11, dif11,"o")
title("Difference between Data and Range for Run 11")
xlabel("wavelength, nm")
ylabel("difference, nm")
figure()
plot(w11,r11)
title("I0 for Run 11")
xlabel("wavelength, nm")
ylabel("I0")
rw11=range(25.0+dw, stop=49.5+dw, length=nel)
lin_interp = LinearInterpolation(rw11, r11)
cub_interp = CubicSplineInterpolation(rw11, r11)
fine11 = collect(range(25.0+dw, stop=49.5+dw, length=10*nel))
linr11 = lin_interp.(fine11)
cubr11 = cub_interp.(fine11)
figure()
plot(fine11,cubr11-linr11)
title("Linear Interpolation Error")
xlabel("wavelength, nm")
ylabel("difference, nm")

# Check error made because of wavelength offset
lin_interp = LinearInterpolation(w11, r11)
roff11 = lin_interp.(rw11)
diff11 = abs.(roff11.-r11)
figure()
semilogy(rw11,diff11)
title("I0 Error due to Wavelength")
xlabel("wavelength, nm")
ylabel("|error|, nm")
