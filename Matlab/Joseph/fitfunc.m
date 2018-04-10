function yhat = fitfunc(b, angle,lam)
% compute expected R given the fit parameters b
%
% For this version, I'll assume the layer thickness is known
% and I'll fit the real and imaginary parts of the index of refraction
si=sindx(lam);
sio2=sio2ndx(lam);
%y2o3=y2o3ndx(lam);
n = [1,b(1)+b(2)*1i,sio2,si];
x = [0,145,18,0];%150 - for ebeam9
sigma = [4,0,0,0];%4 - for ebeam9
fractions=polarfind(lam);
yhat = zeros(size(angle));
for i=1:length(angle)
    yhat(i)=Parratt(n,x,angle(i),fractions,lam,sigma);
end