function WDmax=FindMaximumWassersteinDistanceHMMSequence(BeforeLength,AfterLength,CostFunction,ElementTotal,TBefore,ElementBefore,StateProportionBefore,mode)
% This function generate a series of simulated fitting trajectories with
% same trajectory length as the tested fitting trajectory, and their MWD
% values are calculated and sorted from smallest to largest. This function
% is only used in TestNonequilibrium (1,2,3), in hypothesis testing.

%Input:
%       BeforeLength:Trajectory length of FRET trajectory before condition
%                    change  
%       AfterLength:Trajectory length of FRET trajectory after condition
%                   change 
%       CostFunction:The cost function of Wasserstein distance, a vector
%       ElementTotal:States corresponding to elements in cost function
%       TBefore: Transition matrix of the steady-state system before
%                condition change
%       ElementBefore: FRET efficiency vector of the steady-state system
%                      before condition change
%       StateProportionBefore: The state distribution of FRET trajectories
%                              before condition change
%       mode:mode==1: don't weight the WD value; mode==2: weight the WD
%            value according to the trajectory lengths forming two
%            distributions; mode==3:weight the WD value according to the 
%            trajectory length of rest part of yfitafter.

%Output:
%       WDmax:A series of MWD values, sorted from smallest to largest, a
%             vector

if nargin==7 
    mode=2; %default values:mode=2, enable weighting
end
if length(union(ElementTotal,ElementBefore))~=length(ElementTotal)||BeforeLength==0||AfterLength==0||length(CostFunction)~=length(ElementTotal)||length(ElementBefore)~=length(StateProportionBefore)||length(ElementBefore)~=length(diag(TBefore)) %input parameters check
    WDmax=[];
    return
end
CStateProportionBefore=cumsum(StateProportionBefore); 
CStateProportionBefore=[0,CStateProportionBefore];
n=length(ElementBefore);
Emit=eye(n);
TotalLength=BeforeLength+AfterLength;
WDmax=zeros(1,1000);

for i=1:1000
    R=rand(1); %use the state distribution of FRET trajectories before condition change to randomly determine which state the first data point in simulated fitting trajectory belongs to
    StartStateNumber=find(CStateProportionBefore<=R,1,'last');
    ElementBefore2=ElementBefore;
    TBefore2=TBefore;
    Order=1:n;
    Order(1)=StartStateNumber;
    Order(StartStateNumber)=1;
    Order=[Order(1),sort(Order(2:end))];
    ElementBefore2=ElementBefore2(Order);
    TBefore2=TBefore2(Order,Order);
    [seq,yfit] = hmmgenerate(TotalLength,TBefore2,Emit); %generate simulated fitting trajectory by a given markov model
    for j=1:n
        lo=seq==j;
        yfit(lo)=ElementBefore2(j);
    end

    yfitbefore=yfit(1:BeforeLength);
    yfitafter=yfit(BeforeLength+1:end);
    [~,WD]=FindMaximumWassersteinDistance3(yfitbefore,yfitafter,CostFunction,ElementTotal,mode); %calculate its MWD value
    WDmax(i)=WD;
end
WDmax=sort(WDmax); %sort these MWD values from smallest to largest