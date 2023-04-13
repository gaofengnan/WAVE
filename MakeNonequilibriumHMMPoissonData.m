% MakeNonequilibriumHMMPoissonData
% This script is used to generate simulated FRET trajectory. Those FRET
% trajectories only contain Poisson noise. 
file_pathTest=uigetdir('E:\tirf data','choose the tif name folder that you want to use');
file_pathTest2=[file_pathTest,'\','E'];
mkdir(file_pathTest2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mode controls which kinetic rate coefficients set you're going to use to
% generate data, like follows: 
% mode=1 :non-overlapping 7-state system [A-1]
% mode=2 :partially overlapping 5-state system [B]
% mode=3 :inclusion relationship 4 states system [C]
% mode=4 :inclusion relationship 4 states system [D]
% mode=5 :non-overlapping 4-state system [A-2], with fixed BNEST position,
%         only used in second test 
% mode=6 :inclusion relationship 4 states system [C], with fixed BNEST
%         position, only used in forth test
mode=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1,E1 :HMM parameter set of first steady-state system 
% T2,E2 :HMM parameter set of second steady-state system 
% TranState :states in first steady-state system that can lead to BNEST
%            after condition change
% Rate :the probability of BNEST corresponding to each state in TranState
% FRETLength :the FRET region length
% ChangeFrame :the condition change position
Cyclenumber=100; %the number of FRET trajectories you want to generate
Xd=0.1; %crosstalk coefficient, from donor channel to accepter channel
kabg=5; %accepter channel dark count (background level)
kdbg=5; %accepter channel dark count (background level)
kmol=85.8; %total intensity
CrosstalkLength=150; %the Crosstalk region length
BackgroundLength=150; %the Background region length
if mode==1 %non-overlapping 7-state system [A-1]
    T1=[0.9,0.05,0,0.05;0.05,0.8,0.15,0;0,0.15,0.8,0.05;0.05,0,0.05,0.9];E1=[0.8;0.6;0.4;0.2];
    T2=[0.9,0.05,0.05;0.05,0.9,0.05;0.05,0.05,0.9];E2=[0.9;0.5;0.3];
    Rate=[0.1,0.05]; TranState=[0.2,0.8];
    FRETLength=300; 
    ChangeFrame=150; 
elseif mode==2 %partially overlapping 5-state system [B]
    T1=[0.9,0.05,0,0.05;0.05,0.8,0.15,0;0,0.15,0.8,0.05;0.05,0,0.05,0.9];E1=[0.8;0.6;0.4;0.2];
    T2=[0.9,0.05,0.05;0.1,0.9,0;0.1,0,0.9];E2=[0.9;0.6;0.4];
    Rate=[0.1,0.05]; TranState=[0.2,0.8];
    FRETLength=300;
    ChangeFrame=150;
elseif mode==3 %inclusion relationship 4 states system [C]
    T1=[0.9,0.05,0,0.05;0.05,0.8,0.15,0;0,0.15,0.8,0.05;0.05,0,0.05,0.9];E1=[0.8;0.6;0.4;0.2];
    T2=[0.95,0.05;0.05,0.95];E2=[0.6;0.2];
    Rate=[0.1,0.05]; TranState=[0.2,0.8];
    FRETLength=300;
    ChangeFrame=150;
elseif mode==4 %inclusion relationship 4 states system [D]
    T2=[0.9,0.05,0,0.05;0.05,0.8,0.15,0;0,0.15,0.8,0.05;0.05,0,0.05,0.9];E2=[0.8;0.6;0.4;0.2];
    T1=[0.95,0.05;0.05,0.95];E1=[0.6;0.4];
    Rate=0.05; TranState=0.4;
    FRETLength=300;
    ChangeFrame=150;
elseif mode==5 %non-overlapping 4-state system [A-2], with fixed BNEST position
    T2=[0.97,0.03;0.05,0.95];E2=[0.8;0.4];
    T1=[0.95,0.05;0.05,0.95];E1=[0.6;0.2];
    Rate=1; TranState=[];
    FRETLength=300;
    ChangeFrame=180;
elseif mode==6 %inclusion relationship 4 states system [C], with fixed BNEST position
    T1=[0.9,0.05,0,0.05;0.05,0.8,0.15,0;0,0.15,0.8,0.05;0.05,0,0.05,0.9];E1=[0.8;0.6;0.4;0.2];
    T2=[0.95,0.05;0.05,0.95];E2=[0.6;0.2];
    Rate=1; TranState=[];
    SkipNumber=75; %control the number of trajectories without BNEST
    FRETLength=300;
    ChangeFrame=180;
end

% Rate=0; TranState=[]; %test difference
% wdet=0.065; %difference broadening
% Rate=0.05; TranState=[0.4];
% Rate=0.05; TranState=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Realafterlength=zeros(1,Cyclenumber); %record the length of each trajectory after BNEST, one part of GT
Realbeforestate=zeros(1,Cyclenumber); %record the state that BNEST occur of each trajectory, one part of GT
DataSetInformation=cell(4,Cyclenumber); %record the GT of all trajectories
[mt,nt]=size(TranState);
if isempty(TranState) %if all states in first steady-state system can lead to BNEST
    for testnum=1:Cyclenumber
    if mode==6
        if testnum<=SkipNumber
            Rate=0;
        else
            Rate=1;
        end
    end
%     E3=E1; %construct E difference
%     for i=1:length(E3)
%         di=2*wdet*rand(1)-wdet; %broaden E to E±wdet
%         E3(i)=E3(i)+di;
%     end
    i=0; 
    while rand(1)>Rate&&i<FRETLength-ChangeFrame %generate random BNEST position
        i=i+1;
    end
    if i==0
        i=1; %prevent 0 frame jump
    end
    if ChangeFrame+i>=FRETLength
        FRETRegion1=FRETLength;
        FRETRegion2=0;
    else
        FRETRegion1=ChangeFrame+i;
        FRETRegion2=FRETLength-FRETRegion1;
        disp([testnum,FRETRegion2]);
        Realafterlength(testnum)=FRETRegion2;
    end
    %generate Intensity-time trajectories
    [BinA1,BinD1,state1]=MakeHMMPoissonData(kabg,kdbg,kmol,T1,E1,FRETRegion1,Xd);
    [BinA2,BinD2,state2]=MakeHMMPoissonData(kabg,kdbg,kmol,T2,E2,FRETRegion2,Xd);
    state=[state1,state2];

    BinA=[BinA1,BinA2,random('Poisson',kabg+kmol*Xd,1,CrosstalkLength),random('Poisson',kabg,1,BackgroundLength)];
    BinD=[BinD1,BinD2,random('Poisson',kdbg+kmol,1,CrosstalkLength),random('Poisson',kdbg,1,BackgroundLength)];
    

%     figure(testnum)
%     plot(BinA,'r');
%     hold on;
%     plot(BinD,'b');
%     hold off;
    
    %calculate the FRET efficiency trajectory
    I=1:FRETLength;
    Iprime=(FRETLength+1):(FRETLength+CrosstalkLength);
    bkg=(FRETLength+CrosstalkLength+1):(FRETLength+CrosstalkLength+BackgroundLength);
    [betaA,betaD,IbetaA,IbetaD]=precalculationTimeBin(I,Iprime,bkg,[BinD',BinA']); %calculate the four paramteter that will be used in later calculation
    [J,E]=postcalculationTimeBin([BinD',BinA'],betaA,betaD,IbetaA,IbetaD,I); %calculate the Fisher's information and FRET efficiency
    Savedata=zeros(length(BinA),2); 
    Savedata(I(1):I(end),1)=E; %record calculation results
    Savedata(I(1):I(end),2)=J;
    Realbeforestate(testnum)=state(ChangeFrame+i);
    dE = abs(diff(E));
    dE = sort(dE);
    sigma = dE(round(0.682*length(dE)))/sqrt(2); % calculate the noise level
    if sigma>=0.25 %if the noise level is too high, discard this data
        disp('This data fluctuates too much, discard')
    else
%         figure(Cyclenumber+testnum)
%         plot(E,'g');
%         hold on;
%         plot(state,'b');
%         hold off;
        %save calculation results
        writematrix([BinD',BinA'],[file_pathTest '\' 'test' num2str(testnum) '.txt']);
        writematrix([1,FRETLength,FRETLength+1,FRETLength+CrosstalkLength,FRETLength+CrosstalkLength+1,FRETLength+CrosstalkLength+BackgroundLength],[file_pathTest2 '\' 'test' num2str(testnum) '.txt' ' Region.txt']);
        writematrix(Savedata,[file_pathTest2 '\' 'test' num2str(testnum) '.txt' ' Efficiency.txt']); 
        %record GT
        DataSetInformation{1,testnum}=['test' num2str(testnum) '.txt'];
        DataSetInformation{2,testnum}=Realafterlength(testnum);
        DataSetInformation{3,testnum}=Realbeforestate(testnum);
        DataSetInformation{4,testnum}=state';
    end
    %FRETLength=FRETLength+30; %used to generate trajectories with different lengths
    end
    save([file_pathTest '\' 'DataSetInformation.mat'],'DataSetInformation'); %save GT
elseif mt==1 %if only several states in first steady-state system can lead to BNEST
    for testnum=1:Cyclenumber
    [BinA1,BinD1,state]=MakeHMMPoissonData(kabg,kdbg,kmol,T1,E1,FRETLength,Xd); %first, generate a full length trajectory without BNEST
    if FRETLength>=ChangeFrame+1
    %decompose this full length trajectory
    stateAfter=state(ChangeFrame+1:FRETLength);
    [Z,~,IZ]=unique(stateAfter); 
    Tran=find(diff(IZ)~=0)'; 
    len=length(Tran)+1;
    StateLength=zeros(1,len);
    Ef=zeros(1,len);
    group1=[1,Tran];
    group2=[Tran,length(stateAfter)];
    if ~isempty(Tran)
        if Tran(1)==1
            group1(2)=-1;
        end
    end
    Ef(1)=mean(stateAfter(1)); 
    StateLength(1)=group2(1);
    for i=2:len
        if group1(i)==-1
            Ef(i)=stateAfter(2);
            StateLength(i)=group2(i)-2+1;
        else
            if group1(i)+1>group2(i)
                Ef(i)=stateAfter(group1(i));
                StateLength(i)=group2(i)-group1(i)+1;
            else
                Ef(i)=stateAfter(group1(i)+1); 
                StateLength(i)=group2(i)-group1(i);
            end
        end
    end
    %generate random BNEST position
    replace=0;
    StartLength=0;
    for i=1:len
        if ~isempty(find(TranState==Ef(i), 1))
            lo=find(TranState==Ef(i), 1);
            k=0;
            while rand(1)>Rate(lo)&&k<FRETLength-ChangeFrame
                k=k+1;
            end
            if k==0
                k=1; %prevent 0 frame jump
            end
            RandomLength=k;
            if StateLength(i)>=RandomLength
                replace=1;
                Realbeforestate(testnum)=Ef(i);
                if i==1
                    StartLength=RandomLength;
                else
                    StartLength=sum(StateLength(1:i-1))+RandomLength;
                end
                break
            end
        end
    end
    if ChangeFrame+StartLength+1>FRETLength
        replace=0;
        Realbeforestate(testnum)=0;
    end
    else
    replace=0;
    end
    %replace the FRET trajectory after BNEST position, generate Intensity-time trajectories
    if replace==1
        L=FRETLength-StartLength-ChangeFrame;
        [BinA2,BinD2,state2]=MakeHMMPoissonData(kabg,kdbg,kmol,T2,E2,L,Xd);
        BinA=[BinA1(1:ChangeFrame+StartLength),BinA2,random('Poisson',kabg+kmol*Xd,1,CrosstalkLength),random('Poisson',kabg,1,BackgroundLength)];
        BinD=[BinD1(1:ChangeFrame+StartLength),BinD2,random('Poisson',kdbg+kmol,1,CrosstalkLength),random('Poisson',kdbg,1,BackgroundLength)];
        state1=state(1:ChangeFrame+StartLength);
        state=[state1,state2];
        disp([testnum,length(BinD2)]);
        Realafterlength(testnum)=length(BinD2);
    else
        BinA=[BinA1,random('Poisson',kabg+kmol*Xd,1,CrosstalkLength),random('Poisson',kabg,1,BackgroundLength)];
        BinD=[BinD1,random('Poisson',kdbg+kmol,1,CrosstalkLength),random('Poisson',kdbg,1,BackgroundLength)];
    end
    


%     figure(testnum)
%     plot(BinA,'r');
%     hold on;
%     plot(BinD,'b');
%     hold off;
    %calculate the FRET efficiency trajectory
    I=1:FRETLength;
    Iprime=(FRETLength+1):(FRETLength+CrosstalkLength);
    bkg=(FRETLength+CrosstalkLength+1):(FRETLength+CrosstalkLength+BackgroundLength);
    [betaA,betaD,IbetaA,IbetaD]=precalculationTimeBin(I,Iprime,bkg,[BinD',BinA']); %calculate the four paramteter that will be used in later calculation
    [J,E]=postcalculationTimeBin([BinD',BinA'],betaA,betaD,IbetaA,IbetaD,I); %calculate the Fisher's information and FRET efficiency
    Savedata=zeros(length(BinA),2); 
    Savedata(I(1):I(end),1)=E; %record calculation results
    Savedata(I(1):I(end),2)=J;
    dE = abs(diff(E));
    dE = sort(dE);
    sigma = dE(round(0.682*length(dE)))/sqrt(2); % calculate the noise level
    if sigma>=0.25 %if the noise level is too high, discard this data
        disp('This data fluctuates too much, discard')
    else
%         figure(Cyclenumber+testnum)
%         plot(E,'g');
%         hold on;
%         plot(state,'b');
%         hold off;
        %save calculation results
        writematrix(Savedata,[file_pathTest2 '\' 'test' num2str(testnum) '.txt' ' Efficiency.txt']); %save calculation results
        writematrix([BinD',BinA'],[file_pathTest '\' 'test' num2str(testnum) '.txt']);
        writematrix([1,FRETLength,FRETLength+1,FRETLength+CrosstalkLength,FRETLength+CrosstalkLength+1,FRETLength+CrosstalkLength+BackgroundLength],[file_pathTest2 '\' 'test' num2str(testnum) '.txt' ' Region.txt']);
        %record GT
        DataSetInformation{1,testnum}=['test' num2str(testnum) '.txt'];
        DataSetInformation{2,testnum}=Realafterlength(testnum);
        DataSetInformation{3,testnum}=Realbeforestate(testnum);
        DataSetInformation{4,testnum}=state';
    end
    %FRETLength=FRETLength+30;  %used to generate trajectories with different lengths
    end
    save([file_pathTest '\' 'DataSetInformation.mat'],'DataSetInformation'); %save GT
end