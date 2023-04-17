function T=Tcalculation(C)
% This function is only used in HMM and calculate the transition matrix
% that satisfied the detail balance.
[n,~]=size(C);
Csum=sum(C,2);
X=zeros(n,n);
for i=1:n
    for j=i:n
        X(i,j)=C(i,j)+C(j,i);
        X(j,i)=X(i,j);
    end
end
X2=zeros(n,n);
L=0;
while sum(sum((X2-X).^2))>n^2*10^-8&&L<100
    L=L+1;
    X2=X;
    Xsum=sum(X2,2);
    for i=1:n
        X(i,i)=C(i,i)*(Xsum(i)-X2(i,i))/(Csum(i)-C(i,i));
    end
    for i=1:n
        for j=i+1:n
            a=Csum(i)-C(i,j)+Csum(j)-C(j,i);
            b=Csum(i)*(Xsum(j)-X2(j,i))+Csum(j)*(Xsum(i)-X2(i,j))-(C(i,j)+C(j,i))*(Xsum(j)-X2(j,i)+Xsum(i)-X2(i,j));
            c=-(C(i,j)+C(j,i))*(Xsum(i)-X2(i,j))*(Xsum(j)-X2(j,i));
            X(i,j)=(-b+sqrt(b^2-4*a*c))/(2*a);
            X(j,i)=X(i,j);
        end
    end
end
Xsum=sum(X,2);
T=X./Xsum;
end