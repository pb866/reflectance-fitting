function yhat = fitfuncs(b, angle,lam, si, sio2)
% compute expected R given the fit parameters b
%
% For this version, I'll assume the layer thickness is known
% and I'll fit the real and imaginary parts of the index of refraction
%y2o3=y2o3ndx(lam);
n = [1,b(1)+b(2)*1i,sio2,si];
x = [0,16.04,2.19,0];
sigma = [b(3),0,0,0];
fractions=polarfind(lam);
yhat = zeros(size(angle));
for i=1:length(angle)
    yhat(i)=Parratt(n,x,angle(i),fractions,lam,sigma);
end