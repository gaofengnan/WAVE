function [Ereal,Treal,Pireal]=HMM_Max_Likelihood(Esave,Pisave,Tsave,Asave,Likehoodsave,LikehoodIndexsave,StateNumber)
% This function is only used in HMM, after HMM_Core. It will choose the
% best model between the given HMM model set.

% Input:
%        Esave:A set of FRET effeciency vectors belong to different HMM models
%        Pisave:A set of stationary distributions belong to different HMM models
%        Tsave:A set of transition matrixs belong to different HMM models
%        Asave:A collection that records whether an HMM loop reported an error
%        Likehoodsave:A set of likehood of every HMM loop
%        LikehoodIndexsave:A set of likehood index of every HMM loop
%        StateNumber: The HMM model state number

%Output:
%        Ereal:FRET effeciency vector of best HMM model
%        Treal:Transition matrix of best HMM model
%        Pireal:Stationary distribution vector of best HMM model

Lo=find(Asave==0); %find the locations of HMM loops which  came out an error
Asave(Lo)=[];
if isempty(Asave) %if all HMM loops came out an error
    Ereal=zeros(StateNumber,1); 
    Treal=zeros(StateNumber,StateNumber);
    Pireal=zeros(StateNumber,1);
else %if not all HMM loops came out an error
    Esave(:,Lo)=[]; %use these locations, delete the HMM model in error 
    Pisave(:,Lo)=[];
    Tsave(Lo,:)=[];
    Likehoodsave(:,Lo)=[];
    LikehoodIndexsave(:,Lo)=[];
    Likehoodsave = log(Likehoodsave);
    Like = sum(Likehoodsave+LikehoodIndexsave); %get the real likelihood of HMM model
    [~,Lo]=max(Like); %choose the model has biggest likelihood

    Ereal=Esave(:,Lo);
    Pireal=Pisave(:,Lo);
    Treal=reshape(Tsave(Lo,:),StateNumber,StateNumber);
    [Ereal,I]=sort(Ereal); %reorder the E vector, from smallest to largest
    Treal=Treal(I,I); %reorder the state will also reorder transition matrix and stationary distribution
    Pireal=Pireal(I);
end
end