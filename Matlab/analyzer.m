%analyzer.m
%Make sure that all Sample Runs are of the same wavelength
%Steve Turley and Joseph Muhlestein 2008-2018
close all;
%clear all;

%fill_sample;

directory = '../Feb2018/'; 
%replace / with \ on a PC and visa versa -from a mac
base = 'Feb2018';
suffix = '.log';
filename = [directory base suffix];
global datafile comment date time gain scantype grating filter1 filter2 ...
    mono sample detector samplex sampley samplez
rsum(filename);
exist diode;
if ans==0
[ang diode m3 current] = rdat(directory);
end

% Get data on a specific run
dark8=mean(diode{263});
dark9=mean(diode{264});
dark10=mean(diode{265});

% Compute I0 for this wavelength
% Note problem with run 268. Starts at 29 nm instead of 26!
I0=interp1(ang{96},diode{96},26);

I8=diode{266};
I10=diode{267};
R8=(I8-dark8)/(I0-dark8);
R10=(I10-dark10)/100/(I0-dark8);
A8=ang{266};
A10=ang{267};
semilogy(A8,R8,'r.',A10,R10,'g.')