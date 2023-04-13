function [Esave,Pisave,Tsave,Asave,Likehoodsave,LikehoodIndexsave]=HMM_Core(StateNumber,to,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame,FretEfficiency)
% This function is only used in HMM, and as a core function, calculate a
% set of HMM model by prepared data, StateNumber and time scale(to). It can
% only be used after HMM_Data_Preparation.

% Last update:2023.01.04: add input parameter FretEfficiency to initialize
%                         FRET efficiency

% Input:
%        StateNumber: The HMM model state number
%        to:Initial time scale, the default value is 5
%        txt_num:The number of .txt files in file_path
%        Kmol:The mean total intensity of FRET frame in both channel after
%             correction,calculated in every trajectory's FRET region
%        Kabg:The mean background intensity of Accepter channel, calculated
%             in every trajectory's background region.
%        Kdbg:The mean background intensity of Donor channel, calculated in
%             every trajectory's background region.
%        BinAA:The intensity trace of Accepter channel in FRET region, of
%              every trajectory
%        BinDD:The intensity trace of Donor channel in FRET region, of
%              every trajectory
%        FretFrame:Record the length of FRET region in every trajectory
%        FretEfficiency:Optional parameter, used to assign a value to the
%                       initial FRET efficiency instead of generating it
%                       randomly

%Output:
%        Esave:A set of FRET effeciency vectors belong to different HMM models
%        Pisave:A set of stationary distributions belong to different HMM models
%        Tsave:A set of transition matrixs belong to different HMM models
%        Asave:A collection that records whether an HMM loop reported an error
%        Likehoodsave:A set of likehood of every HMM loop
%        LikehoodIndexsave:A set of likehood index of every HMM loop

% StateNumber=11;
% to=5;
if nargin==10
    if length(FretEfficiency)~=StateNumber
        FretEfficiency=[];
    end
else
    FretEfficiency=[];
end
Esave=zeros(StateNumber,100); %create space for variables
Pisave=zeros(StateNumber,100);
Tsave=zeros(100,StateNumber*StateNumber);
Asave=zeros(1,100);
Likehoodsave=zeros(txt_num,100);
LikehoodIndexsave=zeros(txt_num,100);

% We will try to build 100 HMM models(LL=100),and choose the best model 
% between them. Note that: not every HMM loop will result in an HMM model,
% some loop will result in an error beacuse of bad initial parameters. And
% every HMM model is only a local optimal solution around the initial
% parameters, not a global optimal solution. So it is why we need to
% choose the best model from them, by HMM_Max_Likelihood.
tic;
parfor LL=1:100 %every HMM loop calculation is time consuming, so I use parfor(Parallel computing) here
    
    Me=[]; %create space for variables
    E=[];
    Pi=[];
    TT=[];
%     a=[];
    likehood=zeros(txt_num,1);
    likehoodIndex=zeros(txt_num,1);
    try %use the try statement to prevent program termination if an error occurs in one HMM loop
    [E,T,Pi]=HMM_Algorithm2(StateNumber,to); %get initial parameters
    if ~isempty(FretEfficiency)
        E=FretEfficiency;
    end
    TT=T;
    T=zeros(StateNumber,StateNumber);
    L=0;
    Kmol2=Kmol; %preassign variables to reduce the communication overhead during parallel computing
    Kabg2=Kabg;
    Kdbg2=Kdbg;
    BinAA2=BinAA;
    BinDD2=BinDD;
    FretFrame2=FretFrame;
    while sum(sum((T-TT).^2))>StateNumber^2*10^-8&&L<100 %each HMM loop will have up to 100 iterations, if the number of iterations exceeds 100, or if the parameters converge, the iteration will exit
        L=L+1;
        Ws=zeros(StateNumber,txt_num);
        Es=zeros(StateNumber,txt_num);

        S=E;
        Pai=Pi;
        T=TT;        
        C=zeros(length(S),length(S));
        
        likehood=zeros(txt_num,1); %save likehood of every trajectory in single loop
        likehoodIndex=zeros(txt_num,1); %save likehood index of every trajectory in single loop
        for k=1:txt_num
            kmol=Kmol2(k);
            kabg=Kabg2(k);
            kdbg=Kdbg2(k);
            BinA=BinAA2(k,1:FretFrame2(k));
            BinD=BinDD2(k,1:FretFrame2(k));
            Ptol=zeros(length(S),length(BinA));
            

            for i=1:length(S)
            Et=S(i);
            ka=kmol*Et+kabg;kd=kmol*(1-Et)+kdbg;
            pa=-ka+BinA.*log(ka)-log(2*pi.*BinA)./2-BinA.*log(BinA)+BinA;
            pd=-kd+BinD.*log(kd)-log(2*pi.*BinD)./2-BinD.*log(BinD)+BinD;
            ptol=pa+pd;
            Ptol(i,:)=exp(ptol);
            end
            
            % Calcluate alpha(a) and beta(b) by forward and backward algorithm
            aT=zeros(length(BinA),length(S));
%               P=diag(Ptol(:,1))./sum(Ptol(:,1));
%               aT(1,:)=Pai'*P;
            aT(1,:)=Pai';
            b=zeros(length(S),length(BinA));
            b(:,end)=ones(length(S),1);
            Alength=length(BinA);
            for i=2:length(BinA)
                P=diag(Ptol(:,i))./sum(Ptol(:,i));
                aT(i,:)=aT(i-1,:)*T*P;
                P=diag(Ptol(:,Alength-i+2))./sum(Ptol(:,Alength-i+2));
                b(:,Alength-i+1)=T*P*b(:,Alength-i+2);
            end
            a=aT';
            
%             P=diag([Ptol(1,1),Ptol(2,1)]);
%             b(:,end)=P*Pai;
%             for i=length(BinA)-1:-1:1
%                 P=diag(Ptol(:,i+1))./sum(Ptol(:,i+1));
%                 b(:,i)=T*P*b(:,i+1);
%             end

            % Calculate likehood and update weight vector and effeciency vector
            Total=a.*b;
            Totalsum=sum(Total);
            Total=Total./Totalsum(ones(length(S),1),1);
            likehood(k)=sum(a(:,end));
            likehoodIndex(k)=sum(log(sum(Ptol(:,2:end))));
            Ws(:,k)=sum(Total.*(BinA+BinD-kabg-kdbg),2);
            Es(:,k)=sum(Total.*(BinA-kabg),2)./Ws(:,k);

            % Update count matrix
            for i=1:length(BinA)-1
                aa=a(:,i);
                aa=aa(:,ones(1,length(S)));
                bb=diag(b(:,i+1));
                P=diag(Ptol(:,i+1));
                Ct=(aa.*T)*(bb.*P);
                Ct=Ct./sum(sum(Ct));
                C=C+Ct;
            end
        end
%         Csum=sum(C,2);
%         Csum=Csum(:,ones(1,length(S)));
%         TT=C./Csum;
        
        % From count matrix, we can update transition matrix and stationary
        % distribution vector, the effeciency vector will also be updated
        TT=Tcalculation(C);
        I=eye(length(S));
        A=single(TT'-I);
        Pi=null(A,'r');
        Pi=Pi./sum(Pi);
        E=sum(Es.*Ws,2)./sum(Ws,2);
        
    end
    catch Me
    end
    if isempty(Me) %if an HMM loop do not result in an error, we will record the HMM model it came out
        if isempty(find(ismissing(E)==1, 1))&&length(Pi)==StateNumber %make a second check
        Esave(:,LL)=E; %record the HMM model 
        Pisave(:,LL)=Pi;
        Tsave(LL,:)=TT(:)';
        Asave(LL)=1;
        Likehoodsave(:,LL)=likehood;
        LikehoodIndexsave(:,LL)=likehoodIndex;
        end
    end
end
toc


end