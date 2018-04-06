function f = polarfind(lam)
%POLARFIND take a lambda (in nm) and return the corresponding polarization
%     of the ALS (Beamline 6.3.2)
%     Theoretical data on polarization from the synchrotron and verified by
%     actual polarization measurements, provided by private communication
%     from Eric Gullikson
%
%Authors: Joseph Muhlestein and Steve Turley 2008-18
load PE.txt;
x=PE(:,1);
y=PE(:,2);
xx=1:.1:3;
yy=interp1(log10(x),y,xx,'spline');
xxx=10.^xx;
ev=1239.8./lam;
p=interp1(xxx,yy,ev,'spline');
f=(p+1)/2;