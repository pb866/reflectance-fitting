# process AlF3 data
# need to go to this directory first with als()
using ALS
using Serialization
using PyPlot

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
	errorbar(x, y, vy, fmt="o")
	title("Dark Current for Gain $gain")
	xlabel("run")
	ylabel("mean dark current")
end

dplot(gain7_runs, 7)
dplot(gain8_runs, 8)
dplot(gain9_runs, 9)
dplot(gain10_runs, 10)