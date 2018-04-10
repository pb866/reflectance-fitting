% Check need to have a separate dark current for each run.
load('runs.mat');
%% Gain 7
dr=[759 782 786 846 856];
d=zeros(1,length(dr));
r=zeros(length(dr),1);
sumdr = 0;
sumr = 0;
for run=1:length(dr)
    d(run)=mean(runs(dr(run)).diode);
    r(run)=std(runs(dr(run)).diode);
    sumdr = sumdr + d(run)*r(run);
    sumr = sumr + r(run);
end
figure
errorbar(dr,d,r,'o');
title('Dark Current, Gain 7');
xlabel('Run Number');
ylabel('Signal');
fprintf('Average 7 dark current is %f\n', sumdr/sumr);
dark(7)=sumdr/sumr;
%% Gain 8
dr=[698 704 710 717 731 737 743 750 756 758 763 765 773 784 788 799 ...
    805 815 821 828 834 840 849];
d=zeros(1,length(dr));
r=zeros(length(dr),1);
sumdr = 0;
sumr = 0;
for run=1:length(dr)
    d(run)=mean(runs(dr(run)).diode);
    r(run)=std(runs(dr(run)).diode);
    sumdr = sumdr + d(run)*r(run);
    sumr = sumr + r(run);
end
figure
errorbar(dr,d,r,'o');
title('Dark Current, Gain 8');
xlabel('Run Number');
ylabel('Signal');
fprintf('Average 8 dark current is %f\n', sumdr/sumr);
dark(8)=sumdr/sumr;
%% Gain 9
dr=[700 706 712 719 723 727 733 739 745 752 767 775 793 801 807 ...
    811 817 824 830 836 842 852 860];
d=zeros(1,length(dr));
r=zeros(length(dr),1);
sumdr = 0;
sumr = 0;
for run=1:length(dr)
    d(run)=mean(runs(dr(run)).diode);
    r(run)=std(runs(dr(run)).diode);
    sumdr = sumdr + d(run)*r(run);
    sumr = sumr + r(run);
end
figure
errorbar(dr,d,r,'o');
title('Dark Current, Gain 9');
xlabel('Run Number');
ylabel('Signal');
fprintf('Average gain 9 dark current is %f\n', sumdr/sumr);
dark(9) = sumdr/sumr;
%% Gain 10
dr=[696 702 708 715 721 725 729 735 741 748 754 770 778 780 795 ...
    803 813 819 826 832 838 844 854];
d=zeros(1,length(dr));
r=zeros(length(dr),1);
sumdr = 0;
sumr = 0;
for run=1:length(dr)
    d(run)=mean(runs(dr(run)).diode);
    r(run)=std(runs(dr(run)).diode);
    sumdr = sumdr + d(run)*r(run);
    sumr = sumr + r(run);
end
figure
errorbar(dr,d,r,'o');
title('Dark Current, Gain 10');
xlabel('Run Number');
ylabel('Signal');
fprintf('Average 10 dark current is %f\n', sumdr/sumr);
dark(10)=sumdr/sumr;
%% Save this for later
save('dark.mat','dark');