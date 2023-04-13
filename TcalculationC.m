function [T,RemainStateNumber]=TcalculationC(C)
% This function is only used after CountMatrixCalculation and calculate the
% transition matrix that satisfied the detail balance. Delete states that
% don't obey the detail balance or transition too often. The output
% contains the transition matrix as well as the serial numbers of states
% that have not been deleted.
[n,~]=size(C);
RemainStateNumber=1:n;
Out=0;
Count=0; 
while Out==0&&Count<10
    Count=Count+1;
    Out=1;
    L=0;
    [n,~]=size(C);
    Csum=sum(C,2);
    X=zeros(n,n);
    for i=1:n %initial detail balance matrix
        for j=i:n
            X(i,j)=C(i,j)+C(j,i);
            X(j,i)=X(i,j);
        end
    end
    X2=zeros(n,n);
    while sum(sum((X2-X).^2))>n^2*10^-8&&L<100 %iterated detail balance matrix
        L=L+1;
        X2=X;
        Xsum=sum(X2,2);
        for i=1:n
            X(i,i)=C(i,i)*(Xsum(i)-X2(i,i))/(Csum(i)-C(i,i));
        end
        diX=diag(X);
        lo=find(isinf(diX)|isnan(diX)); %see if any state does not obey the detail balance
        if ~isempty(lo) %if some states don't obey the detail balance, delete them, and start a new round
            RemainStateNumber(lo)=[];
            C(lo,:)=[];
            C(:,lo)=[];
            Out=0;
            break
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
    if Out==1
    Xsum=sum(X,2);
    T=X./Xsum;
    diT=diag(T);
    lo=find(isinf(diT)|isnan(diT)|diT<exp(-1)); %see if any state does not obey the detail balance, or transition too often
    if ~isempty(lo) %if some states don't obey the detail balance or transition too often, delete them, and start a new round
        RemainStateNumber(lo)=[];
        C(lo,:)=[];
        C(:,lo)=[];
        Out=0;
    end
    end
end
if length(T)==1
    T=1;
end
end