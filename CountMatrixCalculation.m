function [State,C]=CountMatrixCalculation(y,Length)
% This function is only used in TestNonequilibrium (1,2,3) and calculate
% the count matrix using fitting trajectories before condition change.

% Input:
%       y: connection of all fitting trajectories before condition change
%       Length: record the lengths of fitting trajectories before condition change
% Output:
%       State: record the states in y
%       C: calculated count matrix 

if sum(Length)~=length(y) %input check
C=[];
return
end
State=sort(unique(y));
[m,n]=size(y);
if m>n
    y=y';
end
C=zeros(length(State),length(State));
for i=1:length(Length) %calculate the count matrix by accumulation
    if Length(i)==0
        continue
    elseif i==1
        ypart=y(1:Length(i));
    else
        ypart=y(sum(Length(1:i-1))+1:sum(Length(1:i)));
    end
    ymatrix=zeros(length(State),length(ypart));
    for j=1:length(State)
        ymatrix(j,:)=ypart==State(j);
    end
    for j=1:Length(i)-1
        C=C+ymatrix(:,j)*(ymatrix(:,j+1)');
    end
end