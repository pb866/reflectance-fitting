function [n, k, sn, sk, lambda] = fitc(log, runs, i, i0, index, dark)
%FIT Fit a series of runs at a given wavelength
%   This version uses automatic rather than explicit weighting
%   Syntax: fit(log, runs, i, i0, index, dark)
%     log LogFile object with information from log file
%     runs array of Run objects with run data for all necessary runs
%     i run(s) for reflected signal
%     i0 run for i0
%     index array of Index objects for Y2O3, Si, SiO2
%     dark array of average dark currents for each gain
rfl = Reflectance(log, runs, i, i0, dark);
rfl = rfl.filter(3,90);
angle = rfl.spect(:,1);
refl = rfl.spect(:,2);
lambda = rfl.lambda;
y2o3Index = index(1).at(lambda);
b0 = [real(y2o3Index), imag(y2o3Index)];
siIndex = index(2).at(lambda);
sio2Index = index(3).at(lambda);
[b, r, J, COVB, mse] = nlinfit(angle, refl, ...
    @(x,angle)fitfunc(x,angle,lambda,siIndex,sio2Index), b0, ...
     'ErrorModel', 'Combined');
figure
% get rid of negative data
refl = refl(refl>0);
angle = angle(refl>0);
yf = fitfunc(b,angle,lambda, siIndex, sio2Index);
semilogy(angle, refl, 'ro', angle, yf, 'b-');
xlabel('Angle (degrees)');
ylabel('reflectance');
title(['Y_2O_3 Run ' num2str(i(1)) ':  lambda:  '...
    num2str(round(lambda,1)) ' nm  n = ' num2str(round(b(1),3))... 
    '  k = ' num2str(round(b(2),5)) ' nm ']);
    legend('data', 'fit');
n=b(1);
k=b(2);
sn=sqrt(COVB(1,1));
sk=sqrt(COVB(2,2));
end