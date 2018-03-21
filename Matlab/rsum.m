function rsum(file)
% RSUM Read in an ALS log file and create arrays with critical data

% Construct file name from directory, base, and suffix
global datafile comment date time gain scantype grating filter1 filter2 ...
    mono sample detector samplex sampley samplez
% Open file and read in the formatted data
%fprintf('..' filesep 'Th0704' filesep 'Th0704.log');
fid = fopen(file);
% TAB = char(009); %ascii code for TAB
format = ['%*s%*s%s%s%s%s%*s%*s%n%*s%*s%*s%*s%*s%*s%*s%s' ...
    '%*s%d%*s%*s%s%s%*s%*s%*s%n%*s%n%*s%*s%*s%*s%*s%n%*s' ...
    '%*n%*s%n%*s%n%*s%n%*s'];
C = textscan(fid, format, 'Delimiter', char(009), 'Headerlines', 2);
fclose(fid);
% Construct useful arrays
datafile = C{1};
comment = C{2};
date = C{3};
time = C{4};
gain = C{5};
scantype = C{6};
grating = C{7};
filter1 = C{8};
filter2 = C{9};
mono = C{10};
sample = C{11};
detector = C{12};
samplex = C{13};
sampley = C{14};
samplez = C{15};