# process AlF3 data
# need to go to this directory first with als()
using ALS
using Serialization

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
