function R=Parratt(n,x,theta,fractions,lambda,sigma)
% calculate the reflectance a multilayer stack
% Written by Steve Turley, 6 July 2006
fractionp=1-fractions; % fraction of p-polarized light
S=sqrt(n.^2-cosd(theta)^2);
k=2*pi/lambda;
C=exp(i*2*S.*x*k);
rs=0;
rp=0;

qz=k*sind(theta); %Debye-Waller rougness correction
eta=exp(-2*qz^2*sigma.^2); %Debye-Waller rougness correction

%qz=k*S; %Nove-Croce rougness correction
for m=length(n)-1:-1:1;
    %eta=exp(-2*qz(m)*qz(m+1)*sigma(m)^2); %Nove-Croce rougness correction
    fs = (S(m)-S(m+1))/(S(m)+S(m+1));
    fp = (n(m+1)^2*S(m)-n(m)^2*S(m+1))/(n(m+1)^2*S(m)+n(m)^2*S(m+1));
    rs=C(m)*(fs*eta(m)+rs*eta(m)^2)/(1+fs*rs*eta(m));
    rp=C(m)*(fp*eta(m)+rp*eta(m)^2)/(1+fp*rp*eta(m));
end
R=(fractionp*abs(rp)^2+fractions*abs(rs)^2);
return;