function [T,Pi,E,Fittrace]=HMM_MainSimple(InitialStateNumber,File_Path,InitialFretEfficiency)
% This function is a simplified version of HMM_Main, also can be called in
% other program. In this function, we use HMM method to explain our FRET
% data, and the HMM validation and refinement part have been removed. All 
% FRET data should be prepared in following format
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
%         -->XXX.txt Region.txt (corresponding to XXX.txt, records FRET
%            region, contains 6 elements, like follows [FRETbegin, FRETend,
%            Crosstalkbegin, Crosstalkend, Backgroundbegin, Backgroundend]) 



%%%%%%%%%%%%%%%%%%%%%%%%%Parameter Input Guide%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
% Input 0 parameter: This function will ask you to choose a folder where
% data is saved, and plot the distribution of this set of data, then ask
% you for a state number(can only be positive integers). The HMM model will
% be constructed by this state number. 

% Only input InitialStateNumber: This function will ask you to choose a
% folder where data is saved, and construct the HMM model by 
% InitialStateNumber. 

% Input InitialStateNumber and File_Path: This function use the File_Path
% you given to prepare data, and construct the HMM model by 
% InitialStateNumber.  

% Input 3 parameters:This function use the InitialFretEfficiency you
% given to initialize FRET efficiency rather than generating it randomly,
% and construct the HMM model by InitialStateNumber and
% InitialFretEfficiency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%Output Parameter  Guide%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T:
%    The transition matrix of HMM model. It is a square matrix whose
%    dimension represent the number of states. Its element T(i,j) represent
%    the probability of the trajectory stay at state i in time t, will
%    transit to state j after one unit of time. 
% Pi:
%    The stationary distribution of HMM model. It can be calculated by
%    transition matrix T, in this equation Pi'*T=Pi'. It represent the
%    proportion of time that the system stays in each state after system 
%    relaxes. In our HMM calculation, we default the system to thermal
%    equilibrium, which means the stationary distribution Pi and the
%    transition matrix T should meet the detailed balance equations:
%    Pi(i)*T(i,j)=Pi(j)*T(j,i)
% E: 
%    The mean FRET efficiency of each state. Because the efficiency of each
%    state is not one value but a gaussian distribution around a mean
%    value. So we use the mean value to represent the FRET efficiency of
%    this state.
% Fittrace:
%          The total fitting trajectory. It is stitched together by fitting
%          trajectories corresponding to all FRET efficiency-time
%          trajectories.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set input parameters
if nargin==3
    file_path=File_Path;
    StateNumber=InitialStateNumber;
    FretEfficiency=InitialFretEfficiency;
elseif nargin==2 %input all 3 parameters
    file_path=File_Path;
    StateNumber=InitialStateNumber;
    FretEfficiency=[];
elseif nargin==1 %only input InitialStateNumber
    file_path=0;
    StateNumber=InitialStateNumber;
    FretEfficiency=[];
else %input 0 parameter
    file_path=0;
    StateNumber=0;
    FretEfficiency=[];
end
    
%% Prepare data
if file_path==0 %if no file_path is input, this function will ask you to choose a folder where data is saved
    file_path=uigetdir('E:\tirf data\','choose the filepath that you want to use'); %choose the folder
end
[StateNumber,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame,Peverypoint,~,~]=HMM_Data_Preparation(StateNumber,file_path); %data preparation function
if StateNumber==0 %error handling 1(StateNumber=0)
    disp('You select a wrong folder, or your input StateNumber=0')
    return
end

%% Construct and HMM model 
to=5; %the default value of to(time scale) is 5 
[Esave,Pisave,Tsave,Asave,Likehoodsave,LikehoodIndexsave]=HMM_Core(StateNumber,to,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame,FretEfficiency); %get a set of HMM model
[Ereal,Treal,Pireal]=HMM_Max_Likelihood(Esave,Pisave,Tsave,Asave,Likehoodsave,LikehoodIndexsave,StateNumber); %choose the best HMM model between the given HMM model set
if sum(Ereal)<10^-8&&sum(sum(Treal))<10^-8&&sum(Pireal)<10^-8 %error handling 2(no HMM model in model set)
    T=[];
    Pi=[];
    E=[];
    Fittrace=[];
    disp('The state number you choose can not pass the first round iteration, please select another one')
    return
end
[Fittrace,~,DistributionDif]=HMM_FitTrace_Get_Distribution(StateNumber,Ereal,Pireal,Treal,Peverypoint,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame);%get fit trace, FRET distribution of every state, and the difference between stationary distribution and the fit trace
E=Ereal;
T=Treal;
Pi=Pireal;
%% Report the model
disp('------------------------------------------------------------------------------')
disp(File_Path)
show1=['Total State Number=',num2str(StateNumber)];
disp(show1)
show2=['FRET Efficiency vector=',num2str(E')];
disp(show2)
show3=['Stationary distribution=',num2str(Pi')];
disp(show3)
show4=['The difference in stationary distribution=',num2str(DistributionDif)];
disp(show4)
save([File_Path '\' 'BestHMMParameters2.mat'],'E','T','Pi'); %save the best model
disp('------------------------------------------------------------------------------')
end