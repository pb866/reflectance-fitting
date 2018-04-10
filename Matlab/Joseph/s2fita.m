function [b, mse] = s2fita(sp)
%S2FITA Fit a reflectance spectrum
%      returns fitted values without uncertainties
%      This version adds an oxide layer on the aluminum
global AlIndex AlF3Index Si3N4Index Al2O3Index
alndx = AlIndex.at(sp.lambda);
alf3ndx = AlF3Index.at(sp.lambda);
si3n4ndx = Si3N4Index.at(sp.lambda);
al2o3ndx = Al2O3Index.at(sp.lambda);
indx=[alf3ndx, alndx, si3n4ndx al2o3ndx];
b0=[real(alf3ndx) imag(alf3ndx) 11 24 2]; % n(AlF3) t(AlF3) t(Al) t(Al2O3)

% Do the fit with a customized fit
wl = sp.spect(:,1);
rf = ones(length(wl),1);
wgt = sp.spect(:,2);
options = statset('MaxIter',1000);
[b, r, J, COVB, mse] = nlinfit(wl,rf, ...
    @(x,angle)fitfunc2(x, angle, sp, indx)./wgt, b0, options);
fprintf('Fitting AlF3 index with Al2O3 layer at lambda=%f\n', ...
    sp.lambda);
for i=1:length(b)
    fprintf('b(%d)=%f +/- %6.2e\n', i, b(i), sqrt(COVB(i,i)));
end
fspect=fitfunc2(b, wl, sp, indx);
figure
semilogy(wl,sp.spect(:,2),'.', wl, fspect, '-');
title(['Sample 2, \lambda = ' num2str(sp.lambda,3) 'nm']);
xlabel('\theta, deg');
ylabel('Reflectance');
end

