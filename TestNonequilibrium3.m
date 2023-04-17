%TestNonequilibrium3
% This script is the main part of WAVE. 

% This method can identify and locate non-equilibrium transitions in FRET
% trajectories that result from the changing conditions. The key
% information we require in advance is the time point of the condition
% change (Changeframe), which is determined by the experimental design. 

% The method's (script's) workflow can be broken down into three steps.
% First, we employ the STaSI approach, along with a hidden Markov model, to
% fit the entire dataset and obtain a discretized fitting trajectory for
% each FRET trajectory. In step two, a maximum Wasserstein distance (MWD)
% analysis dissects the discrete fitting trajectories at potential change
% points with the largest Wasserstein distances to identify the state
% compositions before and after the biomolecule’s non-equilibrium state
% transition (BNEST). In the third and final stage, we use forward- and
% backward algorithms based on the MDL principle to reach a refined
% localization of the non-equilibrium transition position on each FRET 
% trajectory.           

% Note that, All FRET data should be prepared in following format
% XXX folder(Contain all FRET trajectories)
%     -->XXX.txt (records intensity-time trajectory, contains 2 columns 
%        data with same length, the first column is the intensity-time
%        trajectory of donor channel, and the second column is the
%        intensity-time trajectory of accepter channel) 
%     --E folder(Contain all efficiency and region information of all
%               trajectories in FRET folder)
%         -->XXX.txt Efficiency.txt (corresponding to XXX.txt, records FRET
%            efficiency-time trajectory in the first column, with the same
%            length as intensity-time trajectory. In FRET region: 
%            [FRETbegin, FRETend], the recorded FRET efficiency~=0, while 
%            in crosstalk region and background region, the recorded FRET 
%            efficiency=0)
%         -->XXX.txt Region.txt (corresponding to XXX.txt, records region
%            boundaries and contains 6 elements, like follows [FRETbegin,
%            FRETend, Crosstalkbegin, Crosstalkend, Backgroundbegin,
%            Backgroundend]) 

% Output:
%       BNESTAnanysis:a cell table
%       -----First row: record the name of each FRET trajectory (.txt file)
%       -----Second row: record the length of FRET trajectory after BNEST
%       -----Third row: record the last state before BNEST
%       -----Fourth row: record the first state after BNEST
%       -----Fifth row: record the trajectory length after condition change and before BNEST
%       -----Sixth row: record the Length of last segment before BNEST
%       -----Seventh row: record the FRET trajectory after BNEST
%       -----Eighth row: record the fitting trajectory after BNEST

file_path=uigetdir('E:\tirf data','choose the filepath that you saven trace'); %choose the folder that contain Intensity trace, E trace and region information
%file_path=file_pathTest;
file_path2=[file_path,'\','E'];
txt_path_list = dir(strcat(file_path,'\','*.txt')); %read all .txt files under this folder
txt_num=length(txt_path_list); %read the number of files in this folder
y=[];
ME=[];
S=struct([]);
%%%%%%%%%%%%%%%%%%%%%%%%%%Parameter Setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Changeframe: is the most important parameter that need to be given first
%before calculation, its value only depends on the design of the experiment

%MethodSelect: select a method to generate simulation fitting trajectories
%during hypothesis testing. MethodSelect=1, we use HMM method (recommend),
%MethodSelect~=1, we use permutation and cutting method.

%FGratioBackward and FGratioForward: control the ratio of the goodness of
%fit to the cost of model in Backward and Forward algorithm MDL equations,
%respectively.

Changeframe=150; 
MethodSelect=1;
FGratioBackward=2;
FGratioForward=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
TraceLength=zeros(1,txt_num); %TraceLength record the number of data points per data trace
for j = 1:txt_num %read ...Efficiency.txt field one by one
    txt_name = txt_path_list(j).name;% read .txt file name
    try
    data2=load(strcat(file_path2,'\',txt_name,' Efficiency.txt'));
    catch ME
    end
    if ~isempty(ME)
        ME=[];
        try
        data2=load(strcat(file_path2,'\',txt_name(1:length(txt_name)-4),' Efficiency.txt'));
        catch ME
        end
    end
    if ~isempty(ME)
        ME=[];
        try
        data2=load(strcat(file_path2,'\',txt_name(1:length(txt_name)-11),' Denoise Efficiency.txt'));
        catch ME
        end
    end
    if isempty(ME)
        data2(data2(:,1)==0,:)=[]; %remove extra 0s from E trace
        y=[y;data2(:,1)]; %y splices E trace of different data traces
        TraceLength(j)=length(data2(:,1)); %record the number of data points of this E trace
        S(j).ss=data2(:,1); %record this E trace
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%First Step: trajectories fitting%%%%%%%%%%%%%%%%%%%%%%
%% Find state number(STaSI part)
ME=[];
try
    ParameterSet=load([file_path,'\BestHMMParameters2.mat']); %check to see if this data set has already been analyzed
catch ME
end
if ~isempty(ME) %if the data set has not been analyzed
up=max(1,max(y)+0.1);
down=min(0,min(y)-0.1);
[sfit,MDLN]=MDLMulti(S);
MDLN2=sum(MDLN,1);
TotalStateNumber=length(MDLN2);
% k=1:TotalStateNumber;
% k=TotalStateNumber./(TotalStateNumber+k);
% MDLN3=MDLN2.*k;
MDLN3=MDLN2;
[~,lo]=min(MDLN3);
StateNumberChoose=TotalStateNumber-lo+1;
choose=0;
while choose==0
yfit=[];
yfitbefore=[];
yfitafter=[];
AfterLength=zeros(1,txt_num);
BeforeLength=zeros(1,txt_num);
for j=1:txt_num
    Tracefit=sfit(j).ss(:,lo);
    yfit=[yfit;Tracefit];
    if TraceLength(j)<=Changeframe
        yfitbefore=[yfitbefore;Tracefit];
        BeforeLength(j)=TraceLength(j);
    else
        yfitbefore=[yfitbefore;Tracefit(1:Changeframe)];
        yfitafter=[yfitafter;Tracefit(Changeframe+1:end)];
        AfterLength(j)=length(Tracefit(Changeframe+1:end));
        BeforeLength(j)=Changeframe;
    end
end
StartNumber=StateNumberChoose-10;
if StartNumber<=0
    StartNumber=1;
end
EndNumber=StateNumberChoose+10;
if EndNumber>TotalStateNumber
    EndNumber=TotalStateNumber;
end
disp('************************') %show STaSI result
disp(['This trace has ',num2str(StateNumberChoose),' states'])
disp('************************')
figure(1)
MDLNStart=TotalStateNumber-StartNumber+1;
MDLNEnd=TotalStateNumber-EndNumber+1;
plot(StartNumber:EndNumber,MDLN3(MDLNStart:-1:MDLNEnd),'b');
hold on;
plot([StateNumberChoose,StateNumberChoose],[min(MDLN3(MDLNStart:-1:MDLNEnd)),max(MDLN3(MDLNStart:-1:MDLNEnd))],'r');
hold off;
xlim([StartNumber,EndNumber]);
ylim([min(MDLN3(MDLNStart:-1:MDLNEnd)),max(MDLN3(MDLNStart:-1:MDLNEnd))]);
figure(2)
plot(y,'b');
hold on;
plot(yfit,'r');
hold off;
xlim([1,length(yfit)]);
ylim([down,up]);
disp('*********setting*********')
accept=input('Accept this result?(input y) or choose another integer?(input a integer)','s');
disp('*********setting*********')
if accept=='y'
    choose=1;
elseif ~isempty(str2double(accept))&&~isnan(str2double(accept))
    accept=str2double(accept);
    if accept==fix(accept)&&accept>=1&&accept<=TotalStateNumber
        StateNumberChoose=accept;
        lo=TotalStateNumber-StateNumberChoose+1;
    elseif accept<1||accept>TotalStateNumber
        disp('*********error*********')
        disp('input is out of range')
        disp('*********error*********')
    else
        disp('*********error*********')
        disp('input is not an integer')
        disp('*********error*********')
    end
else
    disp('*********error*********')
    disp('input is incorrect')
    disp('*********error*********')
end
end
%% Construct model(HMM part) and fitting trace
choose=0;
while choose==0
disp('*********HMM setting*********')
HMMStateNumber=input('input HMM initial state number(an integer)','s');
disp('*********HMM setting*********')
if isempty(HMMStateNumber)
    disp('*********error*********')
    disp('input is empty')
    disp('*********error*********')
elseif ~isempty(str2num(HMMStateNumber))&&~sum(isnan(str2num(HMMStateNumber)))&&length(HMMStateNumber)==1
    HMMStateNumber=str2num(HMMStateNumber);
    if HMMStateNumber==fix(HMMStateNumber)&&HMMStateNumber>=1&&HMMStateNumber<=TotalStateNumber
        choose=1;
    elseif HMMStateNumber<1||HMMStateNumber>TotalStateNumber
        disp('*********error*********')
        disp('input is out of range')
        disp('*********error*********')
        choose=0;
    else
        disp('*********error*********')
        disp('some input is not an integer')
        disp('*********error*********')
        choose=0;
    end
else
    disp('*********error*********')
    disp('input is incorrect')
    disp('*********error*********')
end
end
FretEfficiency=zeros(HMMStateNumber,1);
lo=TotalStateNumber-HMMStateNumber+1;
yfit=[];
for k=1:txt_num
Tracefit=sfit(k).ss(:,lo);
yfit=[yfit;Tracefit];
end
FretEfficiency=sort(unique(yfit));
[~,~,~,Fittrace]=HMM_MainSimple(HMMStateNumber,file_path,FretEfficiency);
else %if the data set has been analyzed, load the result

E=ParameterSet.E;
T=ParameterSet.T;
Pi=ParameterSet.Pi;
[StateNumber,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame,Peverypoint,~,~]=HMM_Data_Preparation(length(E),file_path);
[Fittrace,~,~]=HMM_FitTrace_Get_Distribution(StateNumber,E,Pi,T,Peverypoint,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame);

end

if isempty(Fittrace)
    disp('HMM fitting error')
    return
end
%% Divide the trajectories into two parts at the condition change position
yfit=Fittrace';
yfitbefore=[];
yfitafter=[];
ybefore=[];
yafter=[];
AfterLength=zeros(1,txt_num);
BeforeLength=zeros(1,txt_num);
for j=1:txt_num
    if j==1
        Tracefit=yfit(1:TraceLength(1));
        Trace=y(1:TraceLength(1));
    else
        Tracefit=yfit(sum(TraceLength(1:j-1))+1:sum(TraceLength(1:j)));
        Trace=y(sum(TraceLength(1:j-1))+1:sum(TraceLength(1:j)));
    end
    if TraceLength(j)<=Changeframe
        yfitbefore=[yfitbefore;Tracefit];
        ybefore=[ybefore;Trace];
        BeforeLength(j)=TraceLength(j);
    else
        yfitbefore=[yfitbefore;Tracefit(1:Changeframe)];
        yfitafter=[yfitafter;Tracefit(Changeframe+1:end)];
        ybefore=[ybefore;Trace(1:Changeframe)];
        yafter=[yafter;Trace(Changeframe+1:end)];
        AfterLength(j)=length(Tracefit(Changeframe+1:end));
        BeforeLength(j)=Changeframe;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%Second Step: MWD analysis%%%%%%%%%%%%%%%%%%%%%%%%
ElementBefore=sort(unique(yfitbefore));
ElementAfter=sort(unique(yfitafter));
StateProportionBefore=zeros(1,length(ElementBefore));
for i=1:length(ElementBefore)
    StateProportionBefore(i)=length(find(yfitbefore==ElementBefore(i)))./length(yfitbefore);
end
StateProportionAfter=zeros(1,length(ElementAfter));
for i=1:length(ElementAfter)
    StateProportionAfter(i)=length(find(yfitafter==ElementAfter(i)))./length(yfitafter);
end
figure(1000)
bar(StateProportionBefore)
xticklabels(num2str(ElementBefore))
figure(2000)
bar(StateProportionAfter)
xticklabels(num2str(ElementAfter))

% MethodSelect==1, use HMM method to optimize state set and determine whether this trajectory has a BNEST
if MethodSelect==1
[ElementBefore,C]=CountMatrixCalculation(yfitbefore,BeforeLength);
[TBefore,RemainStateNumber]=TcalculationC(C);
if length(RemainStateNumber)<length(ElementBefore)
    for i=1:length(ElementBefore)
        if isempty(find(RemainStateNumber==i, 1))
            disp('*******State Distinguish Before Condition Change*******')
            disp(['State:',num2str(ElementBefore(i)),' has been removed in fitting trace before condition change'])
            disp(['The proportion of this state is:',num2str(StateProportionBefore(i))])
            disp('*******************************************************')
        end
    end
    ElementBefore=ElementBefore(RemainStateNumber);
    StateProportionBefore=StateProportionBefore(RemainStateNumber);
    StateProportionBefore=StateProportionBefore./sum(StateProportionBefore);
    figure(1000)
    bar(StateProportionBefore)
    xticklabels(num2str(ElementBefore))
end
end

% Construct Wasserstein cost function
ElementAll=sort(union(ElementAfter,ElementBefore));
StateProportionAfterTotal=zeros(1,length(ElementAll));
StateProportionBeforeTotal=zeros(1,length(ElementAll));
for i=1:length(ElementAll)
    if ~isempty(StateProportionAfter(ElementAfter==ElementAll(i)))
    StateProportionAfterTotal(i)=StateProportionAfter(ElementAfter==ElementAll(i));
    end
    if ~isempty(StateProportionBefore(ElementBefore==ElementAll(i)))
    StateProportionBeforeTotal(i)=StateProportionBefore(ElementBefore==ElementAll(i));
    end
end

CostFunction=StateProportionAfterTotal./StateProportionBeforeTotal;
CostFunction(CostFunction==Inf)=StateProportionAfterTotal(CostFunction==Inf)*length(ybefore);

% Determine whether this trajectory has a BNEST
yWDFitAfter=[];
yWDAfter=[];
LengthWD=zeros(1,txt_num);
for j=1:txt_num
    if BeforeLength(j)==0||AfterLength(j)==0
        yfitbeforepart=[];
        yfitafterpart=[];
        yafterpart=[];
    elseif j==1
        yfitbeforepart=yfitbefore(1:BeforeLength(1));
        yfitafterpart=yfitafter(1:AfterLength(1));
        yafterpart=yafter(1:AfterLength(1));
    else
        yfitbeforepart=yfitbefore(sum(BeforeLength(1:j-1))+1:sum(BeforeLength(1:j)));
        yfitafterpart=yfitafter(sum(AfterLength(1:j-1))+1:sum(AfterLength(1:j)));
        yafterpart=yafter(sum(AfterLength(1:j-1))+1:sum(AfterLength(1:j)));
    end


    if MethodSelect==1 %MethodSelect==1, use HMM method to determine whether this trajectory has a BNEST
        WDmax=FindMaximumWassersteinDistanceHMMSequence(BeforeLength(j),AfterLength(j),CostFunction,ElementAll,TBefore,ElementBefore,StateProportionBefore,2);
        [yfitpart,WD]=FindMaximumWassersteinDistance3(yfitbeforepart,yfitafterpart,CostFunction,ElementAll,2);
        if find(WDmax<WD,1,'last')>=length(WDmax)*(1-0.03)
            yWDFitAfter=[yWDFitAfter;yfitpart];
            yWDAfter=[yWDAfter;yafterpart(end-length(yfitpart)+1:end)];
            LengthWD(j)=length(yfitpart);
        end
    else %MethodSelect~=1, use permutation and cutting method to determine whether this trajectory has a BNEST
        WDmax=FindMaximumWassersteinDistanceOutOrder(BeforeLength(j),AfterLength(j),yfitbefore,CostFunction,ElementAll,2);
        [yfitpart,WD]=FindMaximumWassersteinDistance3(yfitbeforepart,yfitafterpart,CostFunction,ElementAll,2);
        if find(WDmax<WD,1,'last')>=length(WDmax)*(1-0.03)
            yWDFitAfter=[yWDFitAfter;yfitpart];
            yWDAfter=[yWDAfter;yafterpart(end-length(yfitpart)+1:end)];
            LengthWD(j)=length(yfitpart);
        end
    end

    
end

ElementWDAfter=sort(unique(yWDFitAfter));
StateProportionWDAfter=zeros(1,length(ElementWDAfter));
for i=1:length(ElementWDAfter)
    StateProportionWDAfter(i)=length(find(yWDFitAfter==ElementWDAfter(i)))./length(yWDFitAfter);
end
figure(3000)
bar(StateProportionWDAfter)
xticklabels(num2str(ElementWDAfter))

% MethodSelect==1, use HMM method to optimize state set after BNEST
if MethodSelect==1
[ElementWDAfter,C]=CountMatrixCalculation(yWDFitAfter,LengthWD);
RemainStateNumber=1:length(ElementWDAfter);
Tafter=C./sum(C,2); %create transition matrix that does not obey the detail balance
diC=diag(Tafter);
lo=find(isinf(diC)|isnan(diC)|diC<exp(-1));
RemainStateNumber(lo)=[];
if length(RemainStateNumber)<length(ElementWDAfter)
    for i=1:length(ElementWDAfter)
        if isempty(find(RemainStateNumber==i, 1))
            disp('*******State Distinguish After Condition Change*******')
            disp(['State:',num2str(ElementWDAfter(i)),' has been removed in fitting trace after condition change'])
            disp(['The proportion of this state is:',num2str(StateProportionWDAfter(i))])
            disp('*******************************************************')
        end
    end
    ElementWDAfter=ElementWDAfter(RemainStateNumber);
    StateProportionWDAfter=StateProportionWDAfter(RemainStateNumber);
    StateProportionWDAfter=StateProportionWDAfter./sum(StateProportionWDAfter);
    figure(3000)
    bar(StateProportionWDAfter)
    xticklabels(num2str(ElementWDAfter))
end
end

% Warning: Length proportion of trajectory after BNEST is too low
if length(yWDFitAfter)/length(yfitafter)<0.10
    display1='The trace length proportion after find maximum wasserstein distance ';
    display2=['compare to trace length after condition change is:', num2str(length(yWDFitAfter)/length(yfitafter)*100),'%' newline   ...
        'Less than 10%'];
    display3=['It is usually caused by 1)The composition of states before and after ' newline  ...
        'irreversible conformational change is basically the same, and there is ' newline ...
        'no order of magnitude difference in dynamics. 2)The position of ' newline ...
        'irreversible conformational change is too far back in data trace and' newline ...
        'the overall length of data trace after irreversible conformational ' newline...
        'change is too short to be properly recognized'];
    disp(display1)
    disp(display2)
    disp('**************************************************************')
    disp(display3)
    disp('**************************************************************')
    display3=input('Less than 10%, Do you want to continue?(input y to continue)','s');
    disp('**************************************************************')
    if display3~='y'
        disp('Stop program')
        return
    end
end

% Eliminate any states from the state set before or after BNEST if their
% proportion increased by more than 10 times after or before BNEST.
ElementTotal=sort(union(ElementBefore,ElementWDAfter));

ElementAfterUnique=[];
ElementBeforeUnique=[];

for i=1:length(ElementTotal)
    if isempty(find(ElementBefore==ElementTotal(i), 1))
        ElementAfterUnique(end+1)=ElementTotal(i);
    elseif isempty(find(ElementWDAfter==ElementTotal(i), 1))
        ElementBeforeUnique(end+1)=ElementTotal(i);
    elseif StateProportionBefore(ElementBefore==ElementTotal(i))/StateProportionWDAfter(ElementWDAfter==ElementTotal(i))<=1/10
        ElementAfterUnique(end+1)=ElementTotal(i);
    elseif StateProportionBefore(ElementBefore==ElementTotal(i))/StateProportionWDAfter(ElementWDAfter==ElementTotal(i))>=10
        ElementBeforeUnique(end+1)=ElementTotal(i);
    else
        ElementAfterUnique(end+1)=ElementTotal(i);
        ElementBeforeUnique(end+1)=ElementTotal(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%Third Step: Forward/Backward algorithm%%%%%%%%%%%%%%%%%
% Determine the corresponding scenario by comparing state sets before and after BNEST
CommonState=sort(intersect(ElementAfterUnique,ElementBeforeUnique));
AfterLengthResult=zeros(1,txt_num);
yFitResult=[];
% Start forward/backward algorithm
if isempty(CommonState) %Scenario A, completely different state sets appear after BNEST
    LengthBackward=AfterLength-LengthWD;
    LengthForward=LengthWD;
    yfitbackward=[];
    ybackward=[];

    for i=1:length(AfterLength)
        if LengthBackward(i)~=0
            if i==1
                yfitbackward=[yfitbackward;yfitafter(1:LengthBackward(1))];
                ybackward=[ybackward;yafter(1:LengthBackward(1))];
            else
                yfitbackward=[yfitbackward;yfitafter(sum(AfterLength(1:i-1))+1:sum(AfterLength(1:i-1))+LengthBackward(i))];
                ybackward=[ybackward;yafter(sum(AfterLength(1:i-1))+1:sum(AfterLength(1:i-1))+LengthBackward(i))];
            end
        end
    end
    [~,~,Lengthbest]=FindBestTracePart2(yfitbackward,ybackward,LengthBackward,ElementAfterUnique,1,FGratioBackward);
    LengthForward=LengthForward+Lengthbest;
    yfitaforward=[];
    yforward=[];
    for i=1:length(AfterLength)
        if LengthForward(i)~=0
            if i==1
                yfitaforward=[yfitaforward;yfitafter(AfterLength(1)-LengthForward(1)+1:AfterLength(1))];
                yforward=[yforward;yafter(AfterLength(1)-LengthForward(1)+1:AfterLength(1))];
            else
                yfitaforward=[yfitaforward;yfitafter(sum(AfterLength(1:i))-LengthForward(i)+1:sum(AfterLength(1:i)))];
                yforward=[yforward;yafter(sum(AfterLength(1:i))-LengthForward(i)+1:sum(AfterLength(1:i)))];
            end
        end
    end
    [yFitResult,~,AfterLengthResult]=FindBestTracePart2(yfitaforward,yforward,LengthForward,ElementBeforeUnique,2,FGratioForward);

elseif length(CommonState)==length(ElementAfterUnique)&&length(CommonState)==length(ElementBeforeUnique) %Scenario E, state sets before and after BNEST are identical
    disp(['The composition of states before and after irreversible' newline  ...
        'conformational change is basically the same, and there is ' newline ...
        'no order of magnitude difference in dynamics.'])
    AfterLengthResult=LengthWD;
    for i=1:txt_num
        if LengthWD(i)~=0
            AfterLengthResult(i)=LengthWD(i);
            if i==1
                yFitResult=[yFitResult;yWDFitAfter(1:LengthWD(1))];
            else
                yFitResult=[yFitResult;yWDFitAfter(sum(LengthWD(1:i-1))+1:sum(LengthWD(1:i)))];
            end
        end
    end
elseif length(CommonState)==length(ElementAfterUnique)&&length(CommonState)<length(ElementBeforeUnique) %Scenario C, the state set after BNEST is a subset of that before BNEST
    [yfitbest,ybest,Lengthbest]=FindBestTracePart2(yfitafter,yafter,AfterLength,ElementAfterUnique,1,FGratioBackward);
    AfterLengthResult=Lengthbest;
    yFitResult=yfitbest;
elseif length(CommonState)<length(ElementAfterUnique)&&length(CommonState)==length(ElementBeforeUnique) %Scenario D, the state set before BNEST is a subset of that after BNEST
    [yfitbest,ybest,Lengthbest]=FindBestTracePart2(yfitafter,yafter,AfterLength,ElementBeforeUnique,2,FGratioForward);
    AfterLengthResult=Lengthbest;
    yFitResult=yfitbest;
elseif length(CommonState)<length(ElementAfterUnique)&&length(CommonState)<length(ElementBeforeUnique) %Scenario B, two state sets are partially overlapped
    LengthBackward=AfterLength-LengthWD;
    LengthForward=LengthWD;
    yfitbackward=[];
    ybackward=[];

    for i=1:length(AfterLength)
        if LengthBackward(i)~=0
            if i==1
                yfitbackward=[yfitbackward;yfitafter(1:LengthBackward(1))];
                ybackward=[ybackward;yafter(1:LengthBackward(1))];
            else
                yfitbackward=[yfitbackward;yfitafter(sum(AfterLength(1:i-1))+1:sum(AfterLength(1:i-1))+LengthBackward(i))];
                ybackward=[ybackward;yafter(sum(AfterLength(1:i-1))+1:sum(AfterLength(1:i-1))+LengthBackward(i))];
            end
        end
    end
    [~,~,Lengthbest]=FindBestTracePart2(yfitbackward,ybackward,LengthBackward,ElementAfterUnique,1,FGratioBackward);
    LengthForward=LengthForward+Lengthbest;
    yfitaforward=[];
    yforward=[];
    for i=1:length(AfterLength)
        if LengthForward(i)~=0
            if i==1
                yfitaforward=[yfitaforward;yfitafter(AfterLength(1)-LengthForward(1)+1:AfterLength(1))];
                yforward=[yforward;yafter(AfterLength(1)-LengthForward(1)+1:AfterLength(1))];
            else
                yfitaforward=[yfitaforward;yfitafter(sum(AfterLength(1:i))-LengthForward(i)+1:sum(AfterLength(1:i)))];
                yforward=[yforward;yafter(sum(AfterLength(1:i))-LengthForward(i)+1:sum(AfterLength(1:i)))];
            end
        end
    end
    [yFitResult,~,AfterLengthResult]=FindBestTracePart2(yfitaforward,yforward,LengthForward,ElementBeforeUnique,2,FGratioForward);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%Record calculation results and plot figure%%%%%%%%%%%%%%%

if length(AfterLengthResult)==txt_num
    StateBefore=zeros(1,txt_num);
    StateAfter=zeros(1,txt_num);
    TotalTransformLength=-ones(1,txt_num);
    StateTransformLength=-ones(1,txt_num);
%     figure(1)
    if TraceLength(1)<=Changeframe
%         plot(1:TraceLength(1),y(1:TraceLength(1)),'r','LineWidth',1)
%         hold on;
%         plot(1:TraceLength(1),yfit(1:TraceLength(1)),'b','LineWidth',1.5)
%         hold off;
%         xlabel('Frame');ylabel('FRET Efficiency');
    else
%         plot(1:TraceLength(1),y(1:TraceLength(1)),'r','LineWidth',1)
%         hold on;
%         plot(1:TraceLength(1),yfit(1:TraceLength(1)),'b','LineWidth',1.5)
        if AfterLengthResult(1)~=0
            StateBefore(1)=yfit(TraceLength(1)-AfterLengthResult(1));
            StateAfter(1)=yfit(TraceLength(1)-AfterLengthResult(1)+1);
            TotalTransformLength(1)=TraceLength(1)-AfterLengthResult(1)-Changeframe;
            
            Part=yfit(Changeframe:TraceLength(1)-AfterLengthResult(1));
            if length(Part)>=2
                Count=1;
                for i=1:length(Part)-1
                    if Part(end-i)~=Part(end)
                        break
                    else
                        Count=Count+1;
                    end
                end
                StateTransformLength(1)=Count;
            else
                StateTransformLength(1)=1;
            end
            
%             plot(TraceLength(1)-AfterLengthResult(1)+1:TraceLength(1),yfit(TraceLength(1)-AfterLengthResult(1)+1:TraceLength(1)),'g','LineWidth',1.5);
%             hold off;
        else
%             hold off;
        end
    end
%     ylim([0,1]);
    for j=2:txt_num
%         figure(j)
        frame1=sum(TraceLength(1:j-1))+1;
        frame2=sum(TraceLength(1:j));
        frame3=sum(TraceLength(1:j))-AfterLengthResult(j)+1;
        frame4=TraceLength(j)-AfterLengthResult(j)+1;
        if TraceLength(j)<=Changeframe
%             plot(1:TraceLength(j),y(frame1:frame2),'r','LineWidth',1)
%             hold on;
%             plot(1:TraceLength(j),yfit(frame1:frame2),'b','LineWidth',1.5)
%             hold off;
%             xlabel('Frame');ylabel('FRET Efficiency');
        else
%             plot(1:TraceLength(j),y(frame1:frame2),'r','LineWidth',1)
%             hold on;
%             plot(1:TraceLength(j),yfit(frame1:frame2),'b','LineWidth',1.5)
            if AfterLengthResult(j)~=0
                StateBefore(j)=yfit(frame3-1);
                StateAfter(j)=yfit(frame3);
                TotalTransformLength(j)=TraceLength(j)-AfterLengthResult(j)-Changeframe;
                Part=yfit(frame1+Changeframe-1:frame2-AfterLengthResult(j));
                if length(Part)>=2
                    Count=1;
                    for i=1:length(Part)-1
                        if Part(end-i)~=Part(end)
                            break
                        else
                            Count=Count+1;
                        end
                    end
                    StateTransformLength(j)=Count;
                else
                    StateTransformLength(j)=1;
                end
%                 plot(frame4:TraceLength(j),yfit(frame3:frame2),'g','LineWidth',1.5);
                hold off;
            else
                hold off;
            end
        end 
%         ylim([0,1]);
    end
end

if ~isempty(yFitResult)
    BNESTAnanysis=cell(8,txt_num+1);
    BNESTAnanysis{1,1}='Data name';
    BNESTAnanysis{2,1}='Length of FRET trajectory after BNEST';
    BNESTAnanysis{3,1}='State before BNEST';
    BNESTAnanysis{4,1}='State after BNEST';
    BNESTAnanysis{5,1}='Trajectory length after condition change and before BNEST';
    BNESTAnanysis{6,1}='Length of last segment before BNEST';
    BNESTAnanysis{7,1}='FRET trajectory after BNEST';
    BNESTAnanysis{8,1}='Fitting trajectory after BNEST';
    for i=1:txt_num
        BNESTAnanysis{1,i+1}=txt_path_list(i).name;
        BNESTAnanysis{2,i+1}=AfterLengthResult(i);
        BNESTAnanysis{3,i+1}=StateBefore(i);
        BNESTAnanysis{4,i+1}=StateAfter(i);
        BNESTAnanysis{5,i+1}=TotalTransformLength(i);
        BNESTAnanysis{6,i+1}=StateTransformLength(i);
        if TraceLength(i)~=0&&AfterLengthResult(i)~=0
            if i==1
                BNESTAnanysis{7,i+1}=y(TraceLength(i)-AfterLengthResult(i)+1:TraceLength(i));
                BNESTAnanysis{8,i+1}=yFitResult(1:AfterLengthResult(i));
            else
                BNESTAnanysis{7,i+1}=y(sum(TraceLength(1:i))-AfterLengthResult(i)+1:sum(TraceLength(1:i)));
                BNESTAnanysis{8,i+1}=yFitResult(sum(AfterLengthResult(1:i-1))+1:sum(AfterLengthResult(1:i)));
            end
        end
    end
    save([file_path '\' 'BNESTAnanysis.mat'],'BNESTAnanysis'); %save GT
end