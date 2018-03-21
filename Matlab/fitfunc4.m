function yhat = fitfunc4(b, angle, spect, indx)
% compute expected R given the fit parameters b
%
% For this version, I'll assume the layer thickness is known
% and I'll fit the real and imaginary parts of the index of refraction
n = [1,b(1)+b(2)*1i,b(3)+b(4)*1i,indx];
x = [0,12.5,23,0];
sigma = [0.2,0,0];
yhat = zeros(size(angle));
for i=1:length(angle)
    yhat(i)=Parratt(n, x, angle(i), spect.fracs, spect.lambda, sigma);
end