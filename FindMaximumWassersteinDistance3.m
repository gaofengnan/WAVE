function [yfitpart,WD]=FindMaximumWassersteinDistance3(yfitbefore,yfitafter,CostFunction,ElementTotal,mode)
% This function find the maximum Wasserstein distance on given fitting
% trajectory, which contains two parts: yfitbefore and yfitafter, and only
% used in TestNonequilibrium (1,2,3). We move the WD split point from the 
% beginning of yfitafter, and calculate Wasserstein distances between the
% state distributions before and after each split point (The state
% distribution of yfitbefore+part of yfitafter, and the state distribution
% of rest part of yfitafter), and find the MWD split point with maximum
% Wasserstein distance, record the trajectory after this MWD split point as
% well as the corresponding MWD value.

%Input:
%       yfitbefore:Fitting trajectory that corresponding to FRET trajectory
%                  before condition change 
%       yfitafter:Fitting trajectory that corresponding to FRET trajectory
%                 after condition change 
%       CostFunction:The cost function of Wasserstein distance, a vector
%       ElementTotal:States corresponding to elements in cost function
%       mode:mode==1: don't weight the WD value; mode==2: weight the WD
%            value according to the trajectory lengths forming two
%            distributions; mode==3:weight the WD value according to the 
%            trajectory length of rest part of yfitafter.

%Output:
%       yfitpart:The fitting trajectory after MWD split point
%       WD:MWD value

if nargin==4
    mode=2; %default values:mode=2, enable weighting
end
if isempty(yfitafter)||isempty(yfitbefore)||isempty(CostFunction)||isempty(ElementTotal)||length(CostFunction)~=length(ElementTotal) %input parameters check
    yfitpart=[];
    WD=[];
    return
end
[StateLengthB,StateIntensityB]=TraceDecomposition(yfitbefore); %decomposed yfitbefore and yfitafter
[StateLengthA,StateIntensityA]=TraceDecomposition(yfitafter);
n=length(StateIntensityA);
L=length(ElementTotal);
WDs=zeros(1,n);
UB=sort(unique(StateIntensityB));
UA=sort(unique(StateIntensityA));
ElementAll=sort(union(UB,UA));
if length(union(ElementAll,ElementTotal))~=L %ElementTotal should include all the states in yfitbefore and yfitafter
    yfitpart=[];
    WD=[];
    return
end

% calculate the first WD value when WD split point locate at condition
% change position
StateProbabilityB=zeros(1,L);
StateProbabilityA=zeros(1,L);
for i=1:L %calculate state distributions before and after WD split point
    StateProbabilityB(i)=sum(StateLengthB(StateIntensityB==ElementTotal(i)))/sum(StateLengthB);
    StateProbabilityA(i)=sum(StateLengthA(StateIntensityA==ElementTotal(i)))/sum(StateLengthA);
end
SPdifference=StateProbabilityA-StateProbabilityB;
WDs(1)=sum((SPdifference+abs(SPdifference)).*CostFunction)/2; %calculate the first WD value
if mode==2 %weighting
    WDs(1)=WDs(1)*sqrt((sum(StateLengthA)*sum(StateLengthB))/(sum(StateLengthA)+sum(StateLengthB)));
elseif mode==3
    WDs(1)=WDs(1)*sqrt((sum(StateLengthA)*sum(StateLengthA))/(2*sum(StateLengthA)));
end

% calculate the WD value at each candidate WD split point
for j=1:n-1
    StateLengthBC=[StateLengthB,StateLengthA(1:j)]; %yfitbefore+part of yfitafter
    StateIntensityBC=[StateIntensityB,StateIntensityA(1:j)];
    StateLengthAC=StateLengthA(1+j:end); %rest part of yfitafter
    StateIntensityAC=StateIntensityA(1+j:end);
    
    StateProbabilityB=zeros(1,L);
    StateProbabilityA=zeros(1,L);
    for i=1:L %calculate state distributions before and after WD split point
        StateProbabilityB(i)=sum(StateLengthBC(StateIntensityBC==ElementTotal(i)))/sum(StateLengthBC);
        StateProbabilityA(i)=sum(StateLengthAC(StateIntensityAC==ElementTotal(i)))/sum(StateLengthAC);
    end
    SPdifference=StateProbabilityA-StateProbabilityB;
    WDs(j+1)=sum((SPdifference+abs(SPdifference)).*CostFunction)/2; %calculate the WD value
    if mode==2 %weighting
        WDs(j+1)=WDs(j+1)*sqrt((sum(StateLengthAC)*sum(StateLengthBC))/(sum(StateLengthAC)+sum(StateLengthBC)));
    elseif mode==3
        WDs(j+1)=WDs(j+1)*sqrt((sum(StateLengthAC)*(2*sum(StateLengthA)-sum(StateLengthAC)))/(2*sum(StateLengthA)));
    end
end
[WD,lo]=max(WDs); %find the MWD split point with maximum Wasserstein distance
if lo(1)==1 %record the trajectory after this MWD split point
    yfitpart=yfitafter;
else
    yfitpart=yfitafter(sum(StateLengthA(1:lo(1)-1))+1:end); 
end
end

function [StateLength,StateIntensity]=TraceDecomposition(Trace)
% This function is a called function, decomposed Trace into sections
% recorded by StateLength and StateIntensity. They're all one-dimensional
% vectors and elements in same position are correspond. Used in 
% FindMaximumWassersteinDistance3.
[~,~,IZ]=unique(Trace); 
Tran=find(diff(IZ)~=0)'; %get the state translation location
len=length(Tran)+1; %get the number of segment
StateLength=zeros(1,len); %record the length of every segment
StateIntensity=zeros(1,len); %record the state of every segment
group1=[1,Tran];
group2=[Tran,length(Trace)];
if ~isempty(Tran)
    if Tran(1)==1 %if the first translation happened at frame number 1,set it to -1
        group1(2)=-1;
    end
end
StateIntensity(1)=Trace(1); 
StateLength(1)=group2(1);
for i=2:len
    if group1(i)==-1 %if frame number 1 is regarded as a change point,treated it as frame number 2 when record state and state length
        StateIntensity(i)=Trace(2);
        StateLength(i)=group2(i)-2+1;
    else 
        if group1(i)+1>group2(i) %if the end frame is regarded as a change point,treated it as lastframe-1 when record state and state length
            StateIntensity(i)=Trace(group1(i));
            StateLength(i)=group2(i)-group1(i)+1;
        else %normal situation
            StateIntensity(i)=Trace(group1(i)+1); 
            StateLength(i)=group2(i)-group1(i);
        end
    end
end
end