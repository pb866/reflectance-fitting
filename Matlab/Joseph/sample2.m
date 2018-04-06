%sample2.m
%Analyzie data from sample 2
%Steve Turley March 2018
%close all;
clear all;

load('log.mat'); % log file
load('runs.mat'); % all runs

% dark currents
dark(8) = mean(runs(98).diode);
dark(9) = mean(runs(99).diode);
dark(10) = mean(runs(100).diode);
dark(7) = mean(runs(101).diode);

% 18 nm 
s18=Reflectance(log, runs, 119:121, 97, dark);
f18=s18.filter(1.5,60); % keep part that isn't noisy
% figure();

% Try using the data itself as the fitting weights
wgt = f18.spect(:,2);

% Read in the index of refraction data
AlIndex = Index('Al');
Al2O3Index = Index('Al2O3');
AlF3Index = Index('AlF3');
Si3N4Index = Index('Si3N4');

% Be brave and try fitting Al thickness, AlF3 index, and AlF3 thickness
% all at once.
alndx = AlIndex.at(f18.lambda);
alf3ndx = AlF3Index.at(f18.lambda);
si3n4ndx = Si3N4Index.at(f18.lambda);
indx=[alf3ndx, alndx, si3n4ndx];
b0=[real(alf3ndx) imag(alf3ndx) 5.5 33]; % n(AlF3) t(AlF3) t(Al)

% Do the fit with a customized fit @
wl = f18.spect(:,1);
rf = ones(length(wl),1);
[b, r, J, COVB, mse] = nlinfit(wl,rf, ...
    @(x,angle)fitfunc1(x, angle, f18, indx)./wgt, b0);
b
for i=1:length(b)
    fprintf('sigma %d: %6.2e\n', i, sqrt(COVB(i,i)));
end
mse
fspect=fitfunc1(b, wl, f18, indx);
semilogy(wl,f18.spect(:,2),'.', wl, fspect, '-');
title('Sample 2, \lambda = 18 nm');
xlabel('\theta, deg');
ylabel('Reflectance');
