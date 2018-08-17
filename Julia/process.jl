# process AlF3 data
# need to go to this directory first with als()
using ALS
using Serialization
using PyPlot

# If true, the data will be reinitialized
initRead = false
# If true, plot the dark current data
darkPlots = false
# If true, unfit reflectance data is plotted for each run
reflPlots = true

if initRead
  # To initially read everything
  # This will take a long time
  lf = LogFile()
  println("Reading all data files, please be patient...")
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

# wl = wavelength(lf, 100)
# println("wavelength of run 100 is ", wl)

# Compute and compare dark currents
# Gain 7-runs 27, 40, 51, 60, 68, 84, 98, 112, 125, 136, 147, 169, 182, 191
# Gain 8-runs 13, 28, 41, 52, 61, 69, 85, 99, 113, 126, 137, 148, 170, 183, 192
# Gain 9-runs 12, 29, 42, 53, 62, 70, 86, 100, 114, 127, 138, 139, 149, 171, 184, 193
# Gain 10-runs 14, 30, 43, 54, 63, 71, 87, 101, 115, 128, 140, 150, 172, 185, 194
gain7_runs = (27, 40, 51, 60, 68, 84, 98, 112, 125, 136, 147, 169, 182, 191)
gain8_runs = (13, 28, 41, 52, 61, 69, 85, 99, 113, 126, 137, 148, 170, 183, 192)
gain9_runs = (12, 29, 42, 53, 62, 70, 86, 100, 114, 127, 138, 139, 149, 171, 184, 193)
gain10_runs = (14, 30, 43, 54, 63, 71, 87, 101, 115, 128, 140, 150, 172, 185, 194)
mean(run::Int) = sum(alldata[run].det)/length(alldata[run].det)
function std(run::Int)
	avg = mean(run)
	avg2 = sum(abs2, alldata[run].det)/length(alldata[run].det)
	sqrt(avg2-avg^2)
end

function dplot(runList, gain::Int)
	x = runList
	y = mean.(runList)
	vy = std.(runList)
	figure()
	errorbar(x, y, vy, fmt="o", capsize=5)
	title("Dark Current for Gain $gain")
	xlabel("run")
	ylabel("mean dark current")
end

if darkPlots
  dplot(gain7_runs, 7)
  dplot(gain8_runs, 8)
  dplot(gain9_runs, 9)
  dplot(gain10_runs, 10)
end

# Reflectance data

# 41.3 nm
inum = 19
irun = alldata[inum];
i0run = alldata[11];
darkrun = alldata[12];
dark0run = alldata[13];
lambda = wavelength(lf, inum)
r41p3 = Reflectance(irun, i0run, darkrun, dark0run, lambda)
# drop first points
r41p3 = sub(r41p3, min=1.0)

# 45 nm
inum = 20
irun = alldata[inum];
lambda = wavelength(lf, inum)
r45 = Reflectance(irun, i0run, darkrun, dark0run, lambda)
# drop first point
r45 = sub(r45, min=1.0)

# 35 nm
inum = 21
irun = alldata[inum];
lambda = wavelength(lf, inum)
r35 = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=1.0)

# 30 nm
inum = 22
irun = alldata[inum];
lambda = wavelength(lf, inum)
r30 = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=4.0)

# 27 nm
inum = 23
irun = alldata[inum];
lambda = wavelength(lf, inum)
r27 = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=2.0)

# 26 nm
inum = 24
irun = alldata[inum];
lambda = wavelength(lf, inum)
r26 = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=2.0)

# 25 nm
inum = 25
irun = alldata[inum];
lambda = wavelength(lf, inum)
r25 = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=2.0)


# 31 nm with a different order sorter filter
inum = 31
irun = alldata[inum];
i0run = alldata[26];
darkrun = alldata[28];
dark0run = alldata[27];
lambda = wavelength(lf, inum)
last_lambda = i0run.x[length(i0run.x)]
if lambda>last_lambda
    lambda = last_lambda
end
r30b = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=4.0)

# 27 nm with a different order sorter filter
inum = 32
irun = alldata[inum];
lambda = wavelength(lf, inum)
r27b = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=2.0)

# 26 nm with a different order sorter filter
inum = 33
irun = alldata[inum];
lambda = wavelength(lf, inum)
r26b = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=2.0)

# 25 nm with a different order sorter filter
inum = 34
irun = alldata[inum];
lambda = wavelength(lf, inum)
r25b = sub(Reflectance(irun, i0run, darkrun, dark0run, lambda), min=2.0)

# On plot2 below, red does are the first spectrum and green the second
if reflPlots
    plot(r45)
    plot(r41p3)
    plot(r35)
    plot(r30)
    plot(r30)
    plot2(r30, r30b)
    plot(r27)
    plot(r27b)
    plot2(r27,r27b)
    plot(r26)
    plot(r26b)
    plot2(r26,r26b)
    plot(r25)
    plot(r25b)
    plot2(r25,r25b)
end
