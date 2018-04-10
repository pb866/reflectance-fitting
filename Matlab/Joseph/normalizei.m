function [nrefl,top,bottom] = normalizei(lam, datagain, data, datacurrent,...
    I0gain, I0ang, I0, I0current, datadark, I0dark)
%take all the run info and normalize it, interps to get I0
%interp to get the right values for the I0 diode and current readings
I0i = interp1(I0ang, I0, lam);
I0currenti = interp1(I0ang, I0current, lam);

top=((data / 10^datagain) - datadark) ./ datacurrent;
bottom=((I0i / 10^I0gain) - I0dark) ./ I0currenti;

nrefl=top/bottom;