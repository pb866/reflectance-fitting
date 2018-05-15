% Analyze the runs for sample 14

% clear all
close all

% Load the data which was already saved in binary format
% This should have happened in the initialize script
load('log.mat');
load('runs.mat');
load('dark.mat');

y2o3ndx = Index('Y2O3',9);
sindx = Index('Si',8);
sio2ndx = Index('SiO2',8);
index = [y2o3ndx, sindx, sio2ndx];
n=[];
k=[];
sn=[];
sk=[];
th=[];
sth=[];
lambda = [];
t = 0;
%% 3.8 nm
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [697 699 701], 686, index, dark);
%% 3.9 nm
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [703 705 707], 686, index, dark);
%% 4 nm
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [709 711 714], 686, index, dark);
%% 4.1 nm
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [716 718 720], 686, index, dark);
%% 4.2
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [722 724], 686, index, dark);
%% 4.3
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [726 728], 686, index, dark);
%% 4.4
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [730 732 734], 686, index, dark);
%% 5
% Changed to C filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [736 738 740], 755, index, dark);
%% 6.5
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [749 751 753], 755, index, dark);
%% 12
% Changed to Be filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [764 766 768 769], 760, index, dark);
%% 12.5
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [771 772 774 776 777], 760, index, dark);
% t = t+1;
% [n(t), k(t), sn(t), sk(t),  th(t), sth(t), lambda(t)] = ...
%     fitct(log, runs, [771 772 774 776 777], 760, index, dark);
% fprintf('for lambda=%f, th=%f +/- %f\n', lambda(t), th(t), sth(t));
%% 18
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [845 847 848 850 851 853], 781, index, dark);
%% 20
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [855 857 859 861 862], 781, index, dark);
%% 21
% Changed to Al filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [785 787 789 790 791 792 794], 781, index, dark);
% t = t+1;
% [n(t), k(t), sn(t), sk(t),  th(t), sth(t), lambda(t)] = ...
%     fitct(log, runs, [785 787 789 790 791 792 794], 781, index, dark);
% fprintf('for lambda=%f, th=%f +/- %f\n', lambda(t), th(t), sth(t));
%% 22
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [796 797 798 800 802], 781, index, dark);
%% 23
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [804 806], 781, index, dark);
%% 24
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [808 810 812], 781, index, dark);
%% 25
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [814 816 818], 781, index, dark);
%% 25 repeat, bad data. I'm not sure what's going on here.
% t = t+1;
% [n(t), k(t), sn(t), sk(t), lambda(t)] = ...
%     fitc(log, runs, [864 865 866], 781, index, dark);
%% 26
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [820 822 823 825], 781, index, dark);
%% 27
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [827 829 831], 781, index, dark);
%% 28
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [833 835 837], 781, index, dark);
%% 29
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [839 841 843], 781, index, dark);
%% 30
% Mg Filter--not as good of a fit. Use Al
% t = t+1;
% [n(t), k(t), sn(t), sk(t), lambda(t)] = ...
%     fitc(log, runs, [868 870 871 873], 496, index, dark);
%% 30
% Al Filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [875 877 878 880], 781, index, dark);
%% 36
% Mg Filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [886 887 889:905], 496, index, dark);
%% 42
% Mg Filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [909 911], 496, index, dark);
%% 48
% Mg Filter
t = t+1;
[n(t), k(t), sn(t), sk(t), lambda(t)] = ...
    fitc(log, runs, [914 916], 496, index, dark);
%% Summary
figure
cxron=y2o3ndx.at(lambda);
errorbar(lambda,n,sn,'r.');
hold on
plot(lambda, real(cxron), 'b-');
title('Y_2O_3 Sample 14, n');
xlabel('\lambda, nm');
ylabel('n');
legend('fit','CXRO');
figure
errorbar(lambda, k, sk, 'r.');
hold on
plot(lambda, imag(cxron), 'b-');
title('Y_2O_3 Sample 14, k');
xlabel('\lambda, nm');
ylabel('k');
legend('fit','CXRO');