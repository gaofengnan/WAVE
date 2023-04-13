function [yfitbest,ybest,Lengthbest]=FindBestTracePart2(yfitafter,yafter,AfterLength,IncludeState,mode,FGratio)
%This function combine forward and backward algorithm, selected by mode,
%only used in TestNonequilibrium (1,2,3). We use backward algorithm to find
%trace that is most likely to be data trace after BNEST, and use forward
%algorithm to find trace that is most likely to be data trace before BNEST.
%Use parameter:FGratio to adjust the weight of the goodness of fit in MDL
%formula. 

%Input:
%       yfitafter:The connected curve of the fitting traces, a column vector
%       yafter:The connected curve of the original data traces, a column vector
%       AfterLength:Record the length of every original data trace, a vector
%       IncludeState:A state set contains states that should appear in our finding trace
%       mode:If mode==1, we use backward algorithm, and if mode~=1, we use forward algorithm
%       FGratio:Adjust the weight of the goodness of fit in MDL formula.

%Output:
%       yfitbest:The connected curve of our finding parts that are most likely to be data traces after BNEST, a column vector
%       Lengthbest:Record the length of every finding part, a vector

%% Initial setup
if nargin==4 %if there is no sn or mode input
    mode=1; %default values:mode=1, use backward algorithm
    FGratio=2; %default values:FGratio=2 in backward algorithm
elseif nargin==5 %if there is no FGratio input
    if mode==2
        FGratio=4; %default values:FGratio=4 in forward algorithm
    else
        FGratio=2; %default values:FGratio=2 in backward algorithm
    end
end
yfitbest=[]; %initial outputs and intermediate variables
ybest=[];
Lengthbest=zeros(1,length(AfterLength));
TotalElement=sort(unique(yfitafter));
[m,n]=size(IncludeState);
if n>m
    IncludeState=IncludeState';
end
if length(union(TotalElement,IncludeState))==length(IncludeState)&&mode==1 %if there is no unincluded states, which means all the states in TotalElement are included states, we stop this function
    yfitbest=yfitafter;
    ybest=yafter;
    Lengthbest=AfterLength;
    return
elseif length(union(TotalElement,IncludeState))==length(IncludeState)&&mode==2
    yfitbest=[];
    ybest=[];
    Lengthbest=zeros(1,length(AfterLength));
    return
end
%% Find best part that meet our requirement in each original data trace
for j=1:length(AfterLength) %we look for each original data trace and find best part that meet our requirement
    %% Segment classification on original data trace
    if AfterLength(j)==0 %if the length of this data trace is 0, skip this data trace
        continue
    end
    if mode==1 %if mode==1, we use backward algorithm, read this data trace in positive order
        if j==1
            stateAfter=yfitafter(1:AfterLength(1));
            yafterpart=yafter(1:AfterLength(1));
        else
            stateAfter=yfitafter(sum(AfterLength(1:j-1))+1:sum(AfterLength(1:j)));
            yafterpart=yafter(sum(AfterLength(1:j-1))+1:sum(AfterLength(1:j)));
        end
    else %if mode~=1, we use forward algorithm, read this data trace in reverse order
        if j==1
            stateAfter=flipud(yfitafter(1:AfterLength(1)));
            yafterpart=flipud(yafter(1:AfterLength(1)));
        else
            stateAfter=flipud(yfitafter(sum(AfterLength(1:j-1))+1:sum(AfterLength(1:j))));
            yafterpart=flipud(yafter(sum(AfterLength(1:j-1))+1:sum(AfterLength(1:j))));
        end
    end
    d=diff(yafterpart); %calculate the noise level use Haar
    if ~isempty(d)
        d=abs(d);
        d=sort(d);
        sd=d(round(0.682*length(d)))/sqrt(2); %use the white noise distribution and cumulative distribution function to get noise level
    else
        sd=0.01;
    end
    if sd==0 %used for test
        sd=1;
    end
    [~,~,IZ]=unique(stateAfter); 
    Tran=find(diff(IZ)~=0)'; %get the state translation location
    len=length(Tran)+1; %get the number of segment
    StateLength=zeros(1,len); %record the length of every segment
    Ef=zeros(1,len); %record the state of every segment
    group1=[1,Tran];
    group2=[Tran,length(stateAfter)];
    if ~isempty(Tran)
        if Tran(1)==1 %if the first translation happened at frame number 1,set it to -1
            group1(2)=-1;
        end
    end
    Ef(1)=stateAfter(1); 
    StateLength(1)=group2(1);
    for i=2:len
        if group1(i)==-1 %if frame number 1 is regarded as a change point,treated it as frame number 2 when record state and state length
            Ef(i)=stateAfter(2);
            StateLength(i)=group2(i)-2+1;
        else 
            if group1(i)+1>group2(i) %if the end frame is regarded as a change point,treated it as lastframe-1 when record state and state length
                Ef(i)=stateAfter(group1(i));
                StateLength(i)=group2(i)-group1(i)+1;
            else %normal situation
                Ef(i)=stateAfter(group1(i)+1); 
                StateLength(i)=group2(i)-group1(i);
            end
        end
    end

    
    %% Start Forward/Backward algorithm
    SkipStateNumber=0;
    i=1;
    yfitreplace=[];
    while i<=len %include each segment progressively
        if i==len %intercept a small part of FRET trajectory as well as its corresponding fitting trajectory
            stateaftersmall=stateAfter;
            yaftersmall=yafterpart;
        else
            stateaftersmall=stateAfter(sum(StateLength(1:len-i))+1:end);
            yaftersmall=yafterpart(sum(StateLength(1:len-i))+1:end);
        end
        if length(yaftersmall)<=2 %3 frame limit, the included trace should first longer than 3 frame
            SkipStateNumber=SkipStateNumber+1;
            i=i+1;
            continue
        end
        lengthstateafterunique=0;
        for k=1:length(IncludeState)
            lengthstateafterunique=lengthstateafterunique+length(find(abs(stateaftersmall-IncludeState(k))<0.0001));
        end
        lengthstateafterunique=lengthstateafterunique/length(stateaftersmall);
        if lengthstateafterunique<0.5 %if the proportion of data points that don't belong to the included states exceeds 50%, set MDLC=1 to stop the Forward/Backward algorithm
            MDLC=1;
        else
            if ~isempty(yfitreplace)&&i>=4 %yfitreplace is set to prevent excessive computation, when many segments don't belong to the included states
                stateaftersmall(end-sum(StateLength(len-i+4:len))+1:end)=yfitreplace(end-sum(StateLength(len-i+4:len))+1:end);
            end
            [MDLC,yfitother]=MDLCompare(yaftersmall,sd,stateaftersmall,IncludeState,FGratio); %construct a series of fitting trajectories by included states and compare them with original fitting trajectory
        end
        if ~isempty(MDLC) %see if our inclusion state can better describe this small part of FRET trajectory
            [~,lo]=min(MDLC);
            if lo==length(MDLC) %if original fitting trajectory can best describe this small part of FRET trajectory, stop the Forward/Backward algorithm
                break
            else %if constructed fitting trajectory can best describe this small part of FRET trajectory, the Forward/Backward algorithm continues
                yfitreplace=yfitother(:,lo); %update yfitreplace
                i=i+1;
            end
        else %MDLC=[] means original fitting trajectory only contains included states, the Forward/Backward algorithm continues
            yfitreplace=yfitother; %update yfitreplace
            i=i+1; 
        end
    end


    if i==SkipStateNumber+1&&SkipStateNumber~=0 %if the first >3 frame result stop the Forward/Backward algorithm, start the rollback process
        i=SkipStateNumber;
        while i>=1
            stateaftersmall=stateAfter(sum(StateLength(1:len-i))+1:end); %intercept a small part of FRET trajectory as well as its corresponding fitting trajectory
            yaftersmall=yafterpart(sum(StateLength(1:len-i))+1:end);

            lengthstateafterunique=0;
            for k=1:length(IncludeState)
                lengthstateafterunique=lengthstateafterunique+length(find(abs(stateaftersmall-IncludeState(k))<0.0001));
            end
            lengthstateafterunique=lengthstateafterunique/length(stateaftersmall);
            if lengthstateafterunique<0.5 %if the proportion of data points that don't belong to the included states exceeds 50%, set MDLC=1 to stop the Forward/Backward algorithm
                MDLC=1;
            else
                [MDLC,~]=MDLCompare(yaftersmall,sd,stateaftersmall,IncludeState,FGratio); %construct a series of fitting trajectories by included states and compare them with original fitting trajectory
            end

            if ~isempty(MDLC) %see if our inclusion state can better describe this small part of FRET trajectory
                [~,lo]=min(MDLC);
                if lo==length(MDLC) %if original fitting trajectory can best describe this small part of FRET trajectory, the rollback process continues
                    i=i-1;
                    continue
                else %if constructed fitting trajectory can best describe this small part of FRET trajectory, stop the rollback process
                    i=i+1;
                    break
                end
            elseif isempty(MDLC) %MDLC=[] means original fitting trajectory only contains included states, stop the rollback process
                i=i+1;
                break
            end
        end
    end
    if i==0 %we reach the boundary of the FRET trajectory
        i=1;
    end
    %% Record finding part
    Choose=len-i+2;
    if Choose==len+1
        if mode==1 %in backward algorithm, Choose==len+1 means no segment is included in data trace after BNEST
            continue %skip this data trace
        else %in forward algorithm, Choose==len+1 means no segment is included in data trace before BNEST
            yfitbest=[yfitbest;flipud(stateAfter)]; %the total selected original data trace is likely to be data trace after BNEST
            ybest=[ybest;flipud(yafterpart)];
            Lengthbest(j)=AfterLength(j);
            continue
        end
    end
    if isempty(find(IncludeState==Ef(Choose), 1))&&StateLength(Choose)<3
        Choose=Choose+1;
    end
    

    if mode==1 %in backward algorithm, the included segments is our finding part
        if Choose==1 %in backward algorithm, Choose==1 means the total selected original data trace is likely to be data trace after BNEST
            yfitbest=[yfitbest;stateAfter];
            ybest=[ybest;yafterpart];
            Lengthbest(j)=AfterLength(j);
        else %normal condition
            yfitbest=[yfitbest;stateAfter(sum(StateLength(1:Choose-1))+1:end)];
            ybest=[ybest;yafterpart(sum(StateLength(1:Choose-1))+1:end)];
            Lengthbest(j)=sum(StateLength(Choose:end));
        end
    else %in forward algorithm, trace that is on the original data trace except for the included segments is our finding part
        if Choose==1 %in forward algorithm, Choose(lo)==1 means the total selected original data trace is likely to be data trace before BNEST
            continue %skip this data trace
        else %normal condition
            yfitbest=[yfitbest;flipud(stateAfter(1:sum(StateLength(1:Choose-1))))];
            ybest=[ybest;flipud(yafterpart(1:sum(StateLength(1:Choose-1))))];
            Lengthbest(j)=sum(StateLength(1:Choose-1));
        end
    end
    
end