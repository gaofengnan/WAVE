function [Fittrace,PPart,DistributionDif]=HMM_FitTrace_Get_Distribution(StateNumber,Ereal,Pireal,Treal,Peverypoint,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame)
% This function is only used in HMM, calculate the fit trace, the FRET
% distribution of every state and the difference between stationary
% distribution and the fit trace.

%Input:
%        StateNumber: The HMM model state number
%        Ereal:FRET effeciency vector of best HMM model
%        Pireal:Stationary distribution vector of best HMM model
%        Treal:Transition matrix of best HMM model
%        Peverypoint:Record the poisson distribution of every data point in
%                    total trace, the total trace is a combination of all
%                    trajectories
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

%Output:
%        Fittrace:The fit trace of total trace FRET efficiency
%        PPart:The FRET distribution of every state, calculate by combining
%              all data points' poisson distributions belong to this state
%        DistributionDif:The difference between stationary distribution and
%                        the fit trace

% Trace fitting
Fittrace=[];
for k=1:txt_num
    kmol=Kmol(k);
    kabg=Kabg(k);
    kdbg=Kdbg(k);
    BinA=BinAA(k,1:FretFrame(k));
    BinD=BinDD(k,1:FretFrame(k));
    Ptol=zeros(length(Ereal),length(BinA));


    for i=1:length(Ereal)
    Et=Ereal(i);
    ka=kmol*Et+kabg;kd=kmol*(1-Et)+kdbg;
    pa=-ka+BinA.*log(ka)-log(2*pi.*BinA)./2-BinA.*log(BinA)+BinA;
    pd=-kd+BinD.*log(kd)-log(2*pi.*BinD)./2-BinD.*log(BinD)+BinD;
    ptol=pa+pd;
    Ptol(i,:)=exp(ptol);
    end
    
    % Calcluate alpha(a) and beta(b) by forward and backward algorithm
    aT=zeros(length(BinA),length(Ereal));
    aT(1,:)=Pireal';
    for i=2:length(BinA)
        P=diag(Ptol(:,i))./sum(Ptol(:,i));
        aT(i,:)=aT(i-1,:)*Treal*P;
    end
    a=aT';
    b=zeros(length(Ereal),length(BinA));
    b(:,end)=ones(length(Ereal),1);
    % P=diag([Ptol(1,1),Ptol(2,1)]);
    % b(:,end)=P*Pai;
    for i=length(BinA)-1:-1:1
        P=diag(Ptol(:,i+1))./sum(Ptol(:,i+1));
        b(:,i)=Treal*P*b(:,i+1);
    end
    
    % Calcluate the fit trace for every trajectory, and combine them.
    Total=a.*b;
    Totalsum=sum(Total);
    Total=Total./Totalsum(ones(length(Ereal),1),1);
    [~,fittrace]=max(Total);
    fittrace=Ereal(fittrace);
    for i=1:length(fittrace)
        Fittrace(end+1)=fittrace(i);
    end
end

% Get state distribution and calculate difference
PPart=zeros(1101,StateNumber);
Pi3=zeros(StateNumber,1);
for i=1:StateNumber
    Lo=Fittrace==Ereal(i);
    Pi3(i)=length(find(Lo==1))/length(Fittrace);
    PPart(:,i)=sum(Peverypoint(:,Lo),2);
end
% for i=1:StateNumber
%     figure(1)
%     plot(0:0.001:1.1,PPart(:,i).*1000./sum(sum(PPart)),'k');
%     hold on;
% end
DistributionDif=min(sum(abs(Pireal-Pi3).^2)/(StateNumber*10^-4),sum(abs((Pireal-Pi3)./Pireal).^2)/(StateNumber*10^-2)); %DistributionDif is calculated by sum(abs(difference)^2)/(StateNumber*10^-4) or sum(abs(difference/Pi)^2)/(StateNumber*10^-2)
end
