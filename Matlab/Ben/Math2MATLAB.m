% Read in the CSV files created by Mathematica and convert to MATLAB
% binary data.
%
% Steve Turley, 4/3/18
run={};
for i=1:36
    name=['Analysis/run' num2str(i) '.csv'];
    run{i}=csvread(name);
end
wl=load('Analysis/wavelengths.dat');
save('bendata','run','wl');
% For fitting, it looks like the stack Ben used is as follows
% Si (substrate)
% SiO2 (1.9 nm fit)
% Y2O3 (use n and k from CXRO), t=17 nm)
% vacuum
% roughness about 0.7 nm (fit)