function n=sindx(lambdai)
% return index of Si at the given wavelength lambdai in nm
directory='data\';
material='Si';
postfix='.nk';
file=[directory material postfix];
% display(sprintf('opening file %s',file))
[lambda n beta] = textread(file,'%f %f %f',...
    'headerlines',8);
if(lambdai<lambda(1))
    display('WARNING: lambdai<lambda(1) in sindx')
end
% display(sprintf('%d wavelengths',size(lambda,1)))
if(lambdai>lambda(size(lambda,1)))
    display('WARNING lambdai>lambda(MAX) in sindx')
end
nd=interp1(lambda, n, lambdai,'spline');
b=interp1(lambda, beta, lambdai, 'spline');
n=nd+b*1i;