% Analyze the runs for sample 14

% clear all
close all

% Load the data which was already saved in binary format
% This should have happened in Math2MATLAB.m
load('bendata.mat');

y2o3ndx = Index('Y2O3',9);
sindx = Index('Si',8);
sio2ndx = Index('SiO2',8);
index = [y2o3ndx, sindx, sio2ndx];
n=zeros(length(wl),1);
k=zeros(length(wl),1);
sn=zeros(length(wl),1);
sk=zeros(length(wl),1);
th=zeros(length(wl),1);
sth=zeros(length(wl),1);
th1=zeros(length(wl),1);
sth1=zeros(length(wl),1);
sig=zeros(length(wl),1);
ssig=zeros(length(wl),1);
t = 0;
for i=1:length(run)
    [n(i), k(i), sn(i), sk(i), th(i), sth(i), th1(i), sth1(i), ...
        sig(i), ssig(i)] = fita(run{i}, index, wl(i));
end
%% Summary
figure
lam = linspace(min(wl), max(wl),500);
cxron=y2o3ndx.at(lam);
errorbar(wl,n,sn,'r.');
hold on
plot(lam, real(cxron), 'b-');
title('Y_2O_3 Sample 14, n');
xlabel('\lambda, nm');
ylabel('n');
legend('fit','CXRO');
figure
errorbar(wl, k, sk, 'r.');
hold on
plot(lam, imag(cxron), 'b-');
title('Y_2O_3 Sample 14, k');
xlabel('\lambda, nm');
ylabel('k');
legend('fit','CXRO');
figure
errorbar(wl, th, sth, 'r.');
title('Y_2O_3 Sample 14 Thickness');
xlabel('\lambda, nm');
ylabel('thickness, nm');
% gth = th(th>13 & th<20); % good thicknesses
% gsth = sth(th>13 & th<20); % good standard deviations;
% mth=sum(gth./gsth)/sum(1./gsth); % weighted mean
% fprintf('Mean thickness is %f\n',mth);
figure
errorbar(wl, th1, sth1, 'r.');
title('Y_2O_3 Sample 14 SiO_2 Thickness');
xlabel('\lambda, nm');
ylabel('thickness, nm');
figure
errorbar(wl, sig, ssig, 'r.');
title('Y_2O_3 Sample 14 Roughness');
xlabel('\lambda, nm');
ylabel('\sigma, nm');
%% Save data in csv files
csvwrite('n.csv',n);
csvwrite('k.csv',k);