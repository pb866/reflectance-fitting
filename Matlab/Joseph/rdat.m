function runs=rdat(log)
% read all of the data files in log (Log object)
runs(length(log.C{1}))=Run();
for run=1:length(log.C{1})
    runs(run)=Run(log.fullDataFile(run));
end