%combine 80804 and 80805
%Joseph Muhlestein
%10-13-08

% use to do everthing all at once
%analyzer;

%std{1}=stand(fitted);

%test with CXRO
% cxrolam=(1240./data(:,1)).*10;
% cxron=1-data(:,2);
% cxrok=data(:,3);
cxrolam=min(a(4,:)):.1:max(a(4,:));  
cxron=real(y2o3ndx(cxrolam));
cxrok=imag(y2o3ndx(cxrolam));

% stdn=horzcat(stdn4,stdn5);
% stdk=horzcat(stdk4,stdk5);
% for i=1:length(covar)
%     stdn(i)=sqrt(covar{i}(1,1))*10;
%     stdk(i)=sqrt(covar{i}(2,2))*5;
% end
figure
hold on
%errorbar(a(4,:),a(1,:),stdn,'mo')
plot(a(4,:),a(1,:),'go',cxrolam,cxron,'r:')
hold off
title('n vs. lambda')%?
xlabel('lambda (Å)')
ylabel('n')
legend('Measured n', 'CXRO calculated n')
figure
hold on
%errorbar(a(4,:),a(2,:),stdk,'mo')
plot(a(4,:),a(2,:),'bo',cxrolam,cxrok,'r:')
hold off
title('k vs. lambda')
xlabel('lambda (Å)')
ylabel('k')
legend('Measured k', 'CXRO calculated k')