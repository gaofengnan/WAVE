function [MDLC,yfitother]=MDLCompare(yafterpart,sd,yfit,IncludeState,FGratio)
% This function uses the MDL principle, to see if the fitting trajectory
% constructed with the input IncludeState describes yafterpart better than
% yfit does.If yfit only contains states belonging to IncludeState, we
% think IncludeState can best describe yafterpart, otherwise, we replace
% the segments in yfit that belong to states which are not included by
% IncludeState with those of equal lengths but belong to states of
% IncludeState, to construct a series of fitting trajectories {yfit,a}. 
% Next, we calculate the MDL values of {yfit,a} and yfit using MDL formula 
% with the goodness of fit multiplied by a hyper parameter: FGratio. If 
% {yfit,a} has the smallest MDL value, we think IncludeState can best 
% describe yafterpart, otherwise, the yfit can best describe yafterpart. 

% Input:
%       yafterpart: one part of original FRET trajectory
%       sd: standard deviation of noise fluctuation of original FRET
%           trajectory 
%       yfit: the fitting trajectory corresponding to yafterpart
%       IncludeState:A state set contains states that should appear in
%                    fitting trajectory
%       FGratio: hyper parameter, controls the ratio of the goodness of fit
%                to the cost of model in MDL equation
% Output:
%        MDLC: MDL values of every fitting trajectory in yfitother
%        yfitother: a series of fitting trajectories, with the last column
%                   is yfit 

if nargin==4
    FGratio=2; %default values:FGratio=2
end
N = length(yafterpart);
V = max(yafterpart)-min(yafterpart); %get the domain size V
[n,m]=size(IncludeState);
if m>n
    IncludeState=IncludeState';
%     n=m;
end

[C,~,IZ]=unique(yfit); 
if length(union(C,IncludeState))==length(IncludeState) %if yfit only contains states belonging to IncludeState, we think IncludeState can best describe yafterpart
    MDLC=[];
    yfitother=yfit;
    return
end

% decompose yfit
Tran=find(diff(IZ)~=0)'; %get the state translation location
len=length(Tran)+1; %get the number of segment
StateLength=zeros(1,len); %record the length of every segment
Ef=zeros(1,len); %record the state of every segment
group1=[1,Tran];
group2=[Tran,length(yfit)];
if ~isempty(Tran)
    if Tran(1)==1 %if the first translation happened at frame number 1,set it to -1
        group1(2)=-1;
    end
end
Ef(1)=yfit(1); 
StateLength(1)=group2(1);
for i=2:len
    if group1(i)==-1 %if frame number 1 is regarded as a change point,treated it as frame number 2 when record state and state length
        Ef(i)=yfit(2);
        StateLength(i)=group2(i)-2+1;
    else 
        if group1(i)+1>group2(i) %if the end frame is regarded as a change point,treated it as lastframe-1 when record state and state length
            Ef(i)=yfit(group1(i));
            StateLength(i)=group2(i)-group1(i)+1;
        else %normal situation
            Ef(i)=yfit(group1(i)+1); 
            StateLength(i)=group2(i)-group1(i);
        end
    end
end
IncludeStateMatrix=IncludeState(:,ones(1,length(Ef)));
IncludeStateDifference=min(abs(IncludeStateMatrix-Ef));
lo=find(IncludeStateDifference>0.001); %find segments that belong to states which are not included by IncludeState
lo=sort(lo); 
lo2=find(IncludeStateDifference<=0.001); %find segments that belong to states of IncludeState
lo2=sort(lo2);
if isempty(lo) %if yfit only contains states belonging to IncludeState, we think IncludeState can best describe yafterpart
    MDLC=[];
    yfitother=yfit;
    return
end

% if isempty(lo2) 
%     MDLC=1;
%     return
% end

% construct structure G by trend limits and boundary limits
G=struct([]);
Cluster={};
po=1;
ClusterSmall=lo(po);
while po<=length(lo) %cluster adjacent segments that belong to states which are not included by IncludeState
    if po==length(lo)
        Cluster{end+1}=ClusterSmall;
        po=po+1;
    else
        if lo(po)+1==lo(po+1)
            ClusterSmall=[ClusterSmall,lo(po+1)];
            po=po+1;
        else
            Cluster{end+1}=ClusterSmall;
            po=po+1;
            ClusterSmall=lo(po);
        end
    end
end

for i=1:length(Cluster)
    Forward={};
    ClusterSmall=Cluster{i};
    Cpointer=ClusterSmall(1);
    level=1;
    while Cpointer<=ClusterSmall(end)
        if isempty(find(lo==Cpointer-1, 1))
            if isempty(find(lo2==Cpointer-1, 1))
                if isempty(find(lo2==Cpointer+1, 1)) %no limit
                    Forward{level,1}=IncludeState';
                elseif ~isempty(find(lo2==Cpointer+1, 1)) %only right boundary limit
                    if Ef(Cpointer)<=Ef(Cpointer+1)+0.0001
                        Forward{level,1}=IncludeState(IncludeState<=Ef(Cpointer+1)+0.0001)';
                    elseif Ef(Cpointer)>Ef(Cpointer+1)-0.0001
                        Forward{level,1}=IncludeState(IncludeState>=Ef(Cpointer+1)-0.0001)';
                    end
                end
            elseif ~isempty(find(lo2==Cpointer-1, 1))
                if isempty(find(lo2==Cpointer+1, 1)) %only left boundary limit
                    if Ef(Cpointer)<=Ef(Cpointer-1)+0.0001
                        Forward{level,1}=IncludeState(IncludeState<=Ef(Cpointer-1)+0.0001)';
                    elseif Ef(Cpointer)>Ef(Cpointer-1)-0.0001
                        Forward{level,1}=IncludeState(IncludeState>=Ef(Cpointer-1)-0.0001)';
                    end
                elseif ~isempty(find(lo2==Cpointer+1, 1)) %left/right boundary limit
                    if Ef(Cpointer)<=Ef(Cpointer-1)+0.0001 %first, use left boundary limit
                        IS=IncludeState(IncludeState<=Ef(Cpointer-1)+0.0001);
                    elseif Ef(Cpointer)>Ef(Cpointer-1)-0.0001
                        IS=IncludeState(IncludeState>=Ef(Cpointer-1)-0.0001);
                    end
                    if Ef(Cpointer)<=Ef(Cpointer+1)+0.0001 %then, use right boundary limit
                        IS=IS(IS<=Ef(Cpointer+1)+0.0001);
                    elseif Ef(Cpointer)>Ef(Cpointer+1)-0.0001
                        IS=IS(IS>=Ef(Cpointer+1)-0.0001);
                    end
                    Forward{level,1}=IS';
                end
            end
        elseif ~isempty(find(lo==Cpointer-1, 1))
            ForwardLevelTotal=[];
            for j=1:length(Forward(level-1,:))
                ForwardLevelTotal=[ForwardLevelTotal,Forward{level-1,j}];
            end
            if isempty(find(lo2==Cpointer+1, 1)) %only left trend limit
                for j=1:length(ForwardLevelTotal)
                    if Ef(Cpointer)<=Ef(Cpointer-1)+0.0001 %only left trend limit
                        Forward{level,j}=IncludeState(IncludeState<=ForwardLevelTotal(j)+0.0001)';
                    elseif Ef(Cpointer)>Ef(Cpointer-1)-0.0001
                        Forward{level,j}=IncludeState(IncludeState>=ForwardLevelTotal(j)-0.0001)';
                    end
                end
            elseif ~isempty(find(lo2==Cpointer+1, 1)) %left trend limit and right boundary limit
                for j=1:length(ForwardLevelTotal)
                    if Ef(Cpointer)<=Ef(Cpointer-1)+0.0001 %first, use left trend limit
                        IS=IncludeState(IncludeState<=ForwardLevelTotal(j)+0.0001);
                    elseif Ef(Cpointer)>Ef(Cpointer-1)-0.0001
                        IS=IncludeState(IncludeState>=ForwardLevelTotal(j)-0.0001);
                    end
                    if Ef(Cpointer)<=Ef(Cpointer+1)+0.0001 %then, use right boundary limit
                        IS=IS(IS<=Ef(Cpointer+1)+0.0001);
                    elseif Ef(Cpointer)>Ef(Cpointer+1)-0.0001
                        IS=IS(IS>=Ef(Cpointer+1)-0.0001);
                    end
                    Forward{level,j}=IS';
                end
            end
        end
        level=level+1;
        Cpointer=Cpointer+1;
    end
    G(i).gg=Forward;
end

% use structure G to construct a series of fitting trajectories {yfit,a}
for i=1:length(Cluster) 
    Forward=G(i).gg;
    Rode=[];
    for j=1:length(Forward(end,:)) 
        level=length(Forward(:,end));
        RodeSmall=Forward{level,j};
        n=length(RodeSmall);
        Packpos=j;
        while level>1
            level=level-1;
            lengthrecord=[];
            for m=1:length(Forward(level,:))
                if ~isempty(Forward{level,m})
                    lengthrecord=[lengthrecord,length(Forward{level,m})];
                end
            end
            lengthrecord=[0,cumsum(lengthrecord)];
            recordtemp=Packpos-lengthrecord;
            Packpos=find(recordtemp>0,1,'last');
            potemp=recordtemp(Packpos);
            levelrode=Forward{level,Packpos};
            levelrode=levelrode(potemp);
            RodeSmall=[repmat(levelrode,1,n);RodeSmall];
        end
        Rode=[Rode,RodeSmall];
    end
    G(i).gg=Rode';
end

Eftemplate=Ef;
Eftemplate(lo)=0;
for i=1:length(Cluster) 
    ClusterSmall=Cluster{i};
    Rode=G(i).gg;
    [Em,En]=size(Eftemplate);
    [Rn,~]=size(Rode);
    Eftemplate2=zeros(Em*Rn,En);
    for j=1:Em
        Eftemplate2(Rn*(j-1)+1:Rn*j,:)=repmat(Eftemplate(j,:),Rn,1);
        Eftemplate2(Rn*(j-1)+1:Rn*j,ClusterSmall(1):ClusterSmall(end))=Rode;
    end
    Eftemplate=Eftemplate2;
end

[Em,~]=size(Eftemplate);
yfitother=zeros(Em,N);
yfitother(:,1:StateLength(1))=repmat(Eftemplate(:,1),1,StateLength(1));
for i=2:length(StateLength)
    yfitother(:,sum(StateLength(1:i-1))+1:sum(StateLength(1:i)))=repmat(Eftemplate(:,i),1,StateLength(i));
end
yfitother=yfitother';
yfitother=[yfitother,yfit];

% calculate MDL values of yfitother
MDLC=zeros(Em+1,1);
for i=1:Em+1
    [lnDET, ST ,~] = G4th(yfitother(:,i), sd); %get the 4th part of the GL
    F =1/2/sd*sum(abs(yafterpart-yfitother(:,i))); %get the first part of MDL(i),the goodness of the fit
    k=length(unique(yfitother(:,i)));
    GL = k/2*log(1/2/pi)+k*log(V/sd)+ST/2*log(N)+0.5*lnDET; %get the cost of the model
    MDLC(i)=FGratio*F+GL; 
end
    



end

%% Calculate the 4th part of the GL
function [lnDET, ST ,SL] = G4th(sfit, sd)
[Z,~,IZ]=unique(sfit); %Z contains the unique number in sfit from small to large, IZ can be obtained by sfit=Z(IZ). For example, sfit=[1,1,3,3,2,2],Z=[1,2,3],IZ=[1,1,3,3,2,2]
SN=length(Z); %the state number
Tran=find(diff(IZ)~=0); %the state translation location
ST=length(Tran); %the number of state translation
SL=zeros(1,SN); 
CP=zeros(1,ST);
for i=1:SN
   SL(i)=length(IZ(IZ==i)); %the state length
end
for j=1:ST
    CP2=sfit(Tran(j))-sfit(Tran(j)+1); %the difference before and after state translation
    CP(j)=CP2^2;
end
lnDET=sum(log([SL,CP/sd^2])); %calculate the 4th part of the GL
end