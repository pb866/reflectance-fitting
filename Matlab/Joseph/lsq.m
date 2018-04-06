function lsqu = lsq(b, theta, normed, wgt, lam)
% compute expected R given the fit parameters b
%
% For this version, I'll assume the layer thickness is known
% and I'll fit the real and imaginary parts of the index of refraction
yhat = fitfunc(b, theta, lam);
lsqu=sum((yhat-normed).^2./(wgt));