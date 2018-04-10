% Read the data in once and store it
log=LogFile('../Feb09/','Feb09');
for i=1:length(log.C{1})
    runs(i)=Run(log.fullDataFile(i));
end
save('log.mat','log');
save('runs.mat','runs');