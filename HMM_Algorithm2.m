function [E,T,Pi]=HMM_Algorithm2(n,t)
% This function is only used in HMM and calculate initial transform matrix
% T, efficiency vector E and stationary distribution vector Pi, by a given
% state number n and time scale t.
E=zeros(n,1);
T=zeros(n,n);
for i=1:n
    E(i)=randperm(1000,1)/1000;
    T(i,i)=t;
    for j=i+1:n
        T(i,j)=1/(n-1)*randperm(1000,1)/1000;
        T(j,i)=T(i,j);
    end
end
Tsum=sum(T,2);
Tsum=Tsum(:,ones(1,n));
T=T./Tsum;
I=eye(n);
A=single(T'-I);
Pi=null(A,'r');
Pi=Pi./sum(Pi); %Get the stationary distribution Pi, satisfied Pi'*T=Pi'
end
