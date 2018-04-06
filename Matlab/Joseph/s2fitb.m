function [b, mse] = s2fitb(sp)
%S2FITA Fit a reflectance spectrum
%      returns fitted values without uncertainties
%      This version adds an oxide layer on the aluminum
global Al2O3Index AlF3Index Si3N4Index
al2o3ndx = Al2O3Index.at(sp.lambda);
alf3ndx = AlF3Index.at(sp.lambda);
si3n4ndx = Si3N4Index.at(sp.lambda);
indx=si3n4ndx;
b0=[real(alf3ndx) imag(alf3ndx) 11 real(al2o3ndx) imag(al2o3ndx) 24];
% n(AlF3) t(AlF3) n(Al) t(Al)

% Do the fit with a customized fit
wl = sp.spect(:,1);
rf = ones(length(wl),1);
wgt = sp.spect(:,2);
options = statset('MaxIter',1000);
[b, r, J, COVB, mse] = nlinfit(wl,rf, ...
    @(x,angle)fitfunc3(x, angle, sp, indx)./wgt, b0, options);
fprintf('Fitting AlF3 and Al index and thickness at lambda=%f\n', ...
    sp.lambda);
for i=1:length(b)
    fprintf('b(%d)=%f +/- %6.2e\n', i, b(i), sqrt(COVB(i,i)));
end
fspect=fitfunc3(b, wl, sp, indx);
figure
semilogy(wl,sp.spect(:,2),'.', wl, fspect, '-');
title(['Sample 2, \lambda = ' num2str(sp.lambda,3) 'nm']);
xlabel('\theta, deg');
ylabel('Reflectance');
end

