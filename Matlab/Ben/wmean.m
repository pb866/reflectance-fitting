function [avg,sdev] = wmean(vals,dev)
%wmean Calculated weighted mean and adjusted deviation
%   Syntax: [avg, sdev] = wmean(vals, dev)
%       vals: values to be averaged
%       dev: sigma for values to be averaged
%       avg: weighted mean
%       sdev: adjusted standard deviation for this set
avg = sum(vals./dev)/sum(1./dev);
sdev = std(vals)*sqrt(sum(1./dev.^2)/sum(1./dev).^2);
end

