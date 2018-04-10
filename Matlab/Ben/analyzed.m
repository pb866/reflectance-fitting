% Analyze the runs for sample 14, fitting dark current

% clear all
close all

% Load the data which was already saved in binary format
% This should have happened in Math2MATLAB.m
load('bendata.mat');
% To make sure I get all of the data plotted, add 0.1 nm to the wavelength
% if two neighboring runs 
wlast=-1;
for i=1:length(run)
    if wl(i)==wlast
        wl(i) = wl(i-1)+0.01;
    else
        wlast = wl(i);
    end
end
% Runs 33, 26, 27, 18, 14, 11, 10 have bad data
wl(34)=wl(33);
wl(33)=[];
run(33)=[];
wl(26:27)=[];
run(26:27)=[];
wl(23)=[];
run(23)=[];
wl(19)=wl(18); % save first wavelength
wl(18)=[];
run(18)=[];
wl(14)=[];
run(14)=[];
wl(11)=[];
run(11)=[];
wl(10)=[];
run(10)=[];
y2o3ndx = Index('Y2O3',9);
sindx = Index('Si',8);
sio2ndx = Index('SiO2',8);
index = [y2o3ndx, sindx, sio2ndx];
n=zeros(length(wl),1);
k=zeros(length(wl),1);
sn=zeros(length(wl),1);
sk=zeros(length(wl),1);
dark=zeros(length(wl),1);
sdark=zeros(length(wl),1);
for i=1:length(run)
    [n(i), k(i), sn(i), sk(i), dark(i), sdark(i)] = fitd(run{i}, index, wl(i));
end
%% Use average for duplicated 20 nm runs
[n(21),sn(21)]=wmean(n(21:22),sn(21:22));
[k(21),sk(22)]=wmean(k(21:22),sk(21:22));
[dark(21),sdark(21)]=wmean(dark(21:22),sdark(21:22));
n(22)=[];
k(22)=[];
dark(22)=[];
sn(22)=[];
sk(22)=[];
sdark(22)=[];
wl(22)=[];
%% Use average for duplicated 15 nm runs
[n(17),sn(17)]=wmean(n(17:18),sn(17:18));
[k(17),sk(17)]=wmean(k(17:18),sk(17:18));
[dark(17),sdark(17)]=wmean(dark(17:18),sdark(17:18));
n(18)=[];
k(18)=[];
dark(18)=[];
sn(18)=[];
sk(18)=[];
sdark(18)=[];
wl(18)=[];
%% Use average for duplicated 12.5 nm runs
[n(13),sn(13)]=wmean(n(13:14),sn(13:14));
[k(13),sk(13)]=wmean(k(13:14),sk(13:14));
[dark(13),sdark(13)]=wmean(dark(13:14),sdark(13:14));
n(14)=[];
k(14)=[];
dark(14)=[];
sn(14)=[];
sk(14)=[];
sdark(14)=[];
wl(14)=[];
%% Use average for duplicated 5.5 nm runs
[n(2),sn(2)]=wmean(n(2:4),sn(2:4));
[k(2),sk(2)]=wmean(k(2:4),sk(2:4));
[dark(2),sdark(2)]=wmean(dark(2:4),sdark(2:4));
n(3:4)=[];
k(3:4)=[];
dark(3:4)=[];
sn(3:4)=[];
sk(3:4)=[];
sdark(3:4)=[];
wl(3:4)=[];
%% Summary
figure
lam = linspace(min(wl), max(wl), 500);
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
errorbar(wl, dark, sdark, 'r.');
title('Y_2O_3 Sample 14, dark current');
xlabel('\lambda, nm');
ylabel('dark');
%% Save data in csv files
csvwrite('nd.csv',[wl n]);
csvwrite('kd.csv',[wl k]);