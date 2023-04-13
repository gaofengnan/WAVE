function WDmax=FindMaximumWassersteinDistanceOutOrder(BeforeLength,AfterLength,yfitbefore,CostFunction,ElementTotal,mode)
% This function use permutation and cutting method to generate a series of
% simulated fitting trajectories with same trajectory length as the tested
% fitting trajectory, and their MWD values are calculated and sorted from
% smallest to largest. This function is only used in TestNonequilibrium
% (1,2,3), in hypothesis testing.

%Input:
%       BeforeLength:Trajectory length of FRET trajectory before condition
%                    change  
%       AfterLength:Trajectory length of FRET trajectory after condition
%                   change 
%       yfitbefore:The connected curve of the fitting trajectories before
%                  condition change, a column vector.
%       CostFunction:The cost function of Wasserstein distance, a vector
%       ElementTotal:States corresponding to elements in cost function
%       mode:mode==1: don't weight the WD value; mode==2: weight the WD
%            value according to the trajectory lengths forming two
%            distributions; mode==3:weight the WD value according to the 
%            trajectory length of rest part of yfitafter.

%Output:
%       WDmax:A series of MWD values, sorted from smallest to largest, a
%             vector

if nargin==5 
    mode=2; %default values:mode=2, enable weighting
end
if isempty(BeforeLength)||BeforeLength==0||AfterLength==0 %input parameters check
    WDmax=[];
    return
end

[StateLength,StateIntensity]=TraceDecomposition(yfitbefore); %decomposed yfitbefore
L=length(StateLength);
WDmax=zeros(1,1000);

parfor k=1:1000
seq=randperm(L); %rearrange the segments in yfitbefore
StateIntensity2=StateIntensity;
StateLength2=StateLength;
StateIntensity2=StateIntensity2(seq);
StateLength2=StateLength2(seq);

Trace=TraceConstruct(StateLength2,StateIntensity2); %reconstruct the fitting trajectory using the rearranged segments
TraceA=Trace(1:AfterLength);
TraceA=flipud(TraceA);
TraceB=Trace(AfterLength+1:BeforeLength+AfterLength);
TraceB=flipud(TraceB);
[~,WD]=FindMaximumWassersteinDistance3(TraceB,TraceA,CostFunction,ElementTotal,mode); %calculate its MWD value
WDmax(k)=WD;
end
WDmax=sort(WDmax); %sort these MWD values from smallest to largest
end

function [StateLength,StateIntensity]=TraceDecomposition(Trace)
% This function is a called function, decomposed Trace into sections
% recorded by StateLength and StateIntensity. They're all one-dimensional
% vectors and elements in same position are corresponded. Used in 
% FindMaximumWassersteinDistanceOutOrder.
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
function Trace=TraceConstruct(StateLength,StateIntensity)
% This function is a called function, construct Trace by given segment
% information: StateLength and StateIntensity.Used in 
% FindMaximumWassersteinDistanceOutOrder.
Trace=zeros(sum(StateLength),1);
Trace(1:StateLength(1))=StateIntensity(1);
for i=2:length(StateLength)
    Trace(sum(StateLength(1:i-1))+1:sum(StateLength(1:i)))=StateIntensity(i);
end
end