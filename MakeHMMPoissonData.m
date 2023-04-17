function [BinA,BinD,state]=MakeHMMPoissonData(kabg,kdbg,kmol,T,E,Frame,Xd)
% This function is only used to generate trajectory by a given markov
% model. It is a called function in MakeNonequilibriumHMMPoissonData. The 
% intensity of each point on the trajectory follows the Poisson 
% distribution.

% Input:
%        kabg:The mean background intensity of Accepter channel you want to
%             set
%        kdbg:The mean background intensity of Donor channel you want to
%             set
%        kmol:The mean total intensity of FRET frame in both channel you
%             want to set
%        T:Transition matrix of HMM model
%        E:FRET efficiency vector of HMM model
%        Frame:The total number of data point you want to simulate
%        Xd:Crosstalk from channel D to channel A

% Output:
%        BinA:The simulated intensity trajectory of Accepter channel
%        BinD:The simulated intensity trajectory of Donor channel
%        state:The simulated real FRET effeciency
if nargin==6
    Xd=0;
end
if Frame==0
    BinA=[];
    BinD=[];
    state=[];
else
    n=length(E);
    I=eye(n);
    A=T'-I;
    Pi=null(A,'r');
    Pi=Pi./sum(Pi); %Get the stationary distribution Pi, satisfied Pi'*T=Pi'  
    CStateProportion=cumsum(Pi);
    CStateProportion=[0;CStateProportion];
    R=rand(1);
    StartStateNumber=find(CStateProportion<=R,1,'last');
    Order=1:n;
    Order(1)=StartStateNumber;
    Order(StartStateNumber)=1;
    Order=[Order(1),sort(Order(2:end))];
    E=E(Order);
    T=T(Order,Order);
    Emit=eye(n);
    [seq,state] = hmmgenerate(Frame,T,Emit); %generate trajectory by a given markov model
    BinA=zeros(1,Frame);
    BinD=zeros(1,Frame);
    for i=1:n
        lo=find(seq==i);
        BinA(lo)=E(i)*kmol+kabg;
        BinD(lo)=(1-E(i))*kmol+kdbg;
        state(lo)=E(i);
    end
    BinA=real(BinA);
    BinD=real(BinD);
    BinA=BinA+Xd.*(BinD-kdbg);
    parfor i=1:Frame
        BinA(i)=random('Poisson',BinA(i),1,1); %make trajectories follow the Poisson distribution
        BinD(i)=random('Poisson',BinD(i),1,1);
%         SNR=3.5;
%         BinA(i)=normrnd(BinA(i),BinA(i)/SNR/4,1,1); %make trajectories follow the gaussian distribution
%         BinD(i)=normrnd(BinD(i),BinA(i)/SNR/4,1,1);
    end
end
end
