function [StateNumber,txt_num,Kmol,Kabg,Kdbg,BinAA,BinDD,FretFrame,Peverypoint,pall,yp]=HMM_Data_Preparation(StateNumber,file_path)
% This function is only used in HMM, read and prepare data for HMM
% calculation. Your input should contain file_path and StateNumber. If
% StateNumber you input is 0, this function will plot the distribution
% of this group of data and ask you for another StateNumber input.

% Input:
%       StateNumber: The HMM initial state number you set
%       file_path: The folder path you choose ,contain a group of data

% Output:
%        StateNumber: The HMM initial state number you set
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
%        Peverypoint:Record the poisson distribution of every data point in
%                    total trace, the total trace is a combination of all
%                    trajectories
%        pall:The FRET distribution of this group of data, calculate by
%             combining all data points' poisson distributions
%        yp: Record the FRET efficiency of the total trace

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Note%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The file_path you choose should be the folder made by readdata2. If the
% folder has no .txt file, the output StateNumber will be set to 0.

% In this function, we parpare data by using Poisson distribution. This
% method is essentially a recalculation of the FRET efficiency from
% Intensity-frame trajectories. If the quality of the trajectories is high,
% the result of the recalculated FRET efficiency will be substantially
% consistent with the original result, but if the quality of the track is
% low, the corresponding recalculated FRET efficiency will be significantly
% different compared with original result. In this case, we will look for
% the difference between original result and recalculated result for each
% trace, and if the difference is significant, we will construct a virtual
% intensity-frame trajectory using the original result, and then recalculate
% FRET efficiency from that trajectory.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file_path=uigetdir('E:\tirf data\','choose the filepath that you want to use'); %choose the folder
txt_path_list = dir(strcat(file_path,'\','*.txt')); %read all .txt files under this folder
file_path2=strcat(file_path,'\','E'); %file_path contains a E folder, the path of this E folder is assigned to file_path2
txt_num=length(txt_path_list); %txt_num is the number of .txt files in file_path
if txt_num > 0
    dif=zeros(txt_num,1);
    binA=zeros(txt_num,1);
    binD=zeros(txt_num,1);
    RegionData=zeros(txt_num,6);
    y=[]; %record original FRET efficiency data
    for j = 1:txt_num %read .txt file one by one
            txt_name = txt_path_list(j).name;%read .txt file name
            data=load(strcat(file_path,'\',txt_name));%load .txt file
            [m1,~]=size(data);
            [~,n2]=size(binA);
            if m1>n2
                binA=[binA,zeros(txt_num,m1-n2)];
                binD=[binD,zeros(txt_num,m1-n2)];
            end
            binD(j,1:m1)=data(:,1)'; %D-> Donor Intensity; A-> Acceptor Intensity 
            binA(j,1:m1)=data(:,2)'; 
            ME=[];
            
            % Load region data(data3) and original FRET efficiency data(data2) of every trajectory(data)
            try
            data2=load(strcat(file_path2,'\',txt_name,' Efficiency.txt'));
            data3=load(strcat(file_path2,'\',txt_name,' Region.txt')); %data3 contains the region data of every trajectory
            catch ME
            end
            if ~isempty(ME)
                ME=[];
                try
                data2=load(strcat(file_path2,'\',txt_name(1:length(txt_name)-4),' Efficiency.txt'));
                data3=load(strcat(file_path2,'\',txt_name(1:length(txt_name)-4),' Region.txt'));
                catch ME
                end
            end
            if ~isempty(ME)
                ME=[];
                try
                data2=load(strcat(file_path2,'\',txt_name(1:length(txt_name)-11),' Denoise Efficiency.txt'));
                data3=load(strcat(file_path2,'\',txt_name(1:length(txt_name)-11),' Denoise Region.txt'));
                catch ME
                end
            end
            if isempty(ME)
                y=[y;data2(data3(1):data3(2),1)];
                RegionData(j,:)=data3;
            end
    end
    FretFrame=RegionData(:,2)-RegionData(:,1)+1; %record the length of FRET region in every trajectory
    Fmax=max(FretFrame);
    BinAA=zeros(txt_num,Fmax);
    BinDD=zeros(txt_num,Fmax);
    Kmol=zeros(txt_num,1);
    Kabg=zeros(txt_num,1);
    Kdbg=zeros(txt_num,1);
    Peverypoint=zeros(1101,sum(FretFrame)); %record the distribution of every data point 
    pall=zeros(1101,1); %record the total distribution
    yp=[]; %record recalculated FRET efficiency data
    FretFramePosition=zeros(txt_num,2); %record the starting and ending positions of each trace in total trajectory
    FretFramePosition(1,:)=[1,FretFrame(1)];
    for j=2:txt_num
        FretFramePosition(j,:)=[sum(FretFrame(1:j-1))+1,sum(FretFrame(1:j))]; %record the position of every trajectory in total trace
    end
    for j = 1:txt_num
        data2=y(FretFramePosition(j,1):FretFramePosition(j,2));
        data3=RegionData(j,:);
        regionB=data3(5:6);
        regionC=data3(3:4);
        regionF=data3(1:2);
        bkga=sum(binA(j,regionB(1):regionB(end)))/(regionB(end)-regionB(1)+1); %calculate the background of each channel
        bkgd=sum(binD(j,regionB(1):regionB(end)))/(regionB(end)-regionB(1)+1);
        Ia=sum(binA(j,regionF(1):regionF(end)))/(regionF(end)-regionF(1)+1); %calculate the intensity in FRET region of each channel
        Id=sum(binD(j,regionF(1):regionF(end)))/(regionF(end)-regionF(1)+1);
        IbetaD=sum(binD(j,regionC(1):regionC(end)))/(regionC(end)-regionC(1)+1); %calculate IbetaD
        Ica=sum(binA(j,regionC(1):regionC(end)))/(regionC(end)-regionC(1)+1);%calculate the intensity in CrossTalk region of accepter channel
        Xd=(Ica-bkga)/(IbetaD-bkgd); %calculate the parameter of CrossTalk
        if Xd<0 %if CrossTalk parameter is less than 0, alert user and set it to 0
            Xd=0;
        end
        P=(IbetaD-Id)/(Ia-Ica); %calculate intermediate parameter P
        Ia0=((IbetaD-bkgd)+Xd*(IbetaD-bkgd)*P)/P; %calculate the acctpter real intensity when E=1
        IbetaA=Ia0+bkga; %calculate IbetaA
        betaAD=(IbetaA-bkga)/(IbetaD-bkgd);
        BinA=binA(j,data3(1):data3(2));
        BinD=binD(j,data3(1):data3(2));
        kabg=bkga;
        kdbg=bkgd;
        kmol=mean((BinA-bkga-Xd.*(BinD-bkgd))./betaAD+(BinD-bkgd)); %calculate the mean total intensity of FRET frame in both channel after correction
        BinA=(BinA-bkga-Xd.*(BinD-bkgd))./betaAD+bkga; %correct the trajectory in Accepter channel
        if kmol*-0.1+kdbg<1 %for some trajectories with big SNR, the situation kd=kmol*(1-E)+kdbg<0 will occur, which is not allowed. We add the background level to prevent this from happening  
            Add=-(kmol*-0.1+kdbg)+1;
            kdbg=kdbg+Add;
            BinD=BinD+Add;
        end
        if ~isempty(find(BinA<=0, 1)) %BinA<=0 is not allowed, we add the background level to prevent this from happening
            Add=(-min(BinA))+1;
            kabg=kabg+Add;
            BinA=BinA+Add;
        end
        if ~isempty(find(BinD<=0, 1)) %BinD<=0 is not allowed, we add the background level to prevent this from happening
            Add=(-min(BinD))+1;
            kdbg=kdbg+Add;
            BinD=BinD+Add;
        end
        ps=zeros(length(BinA),1101);
        p=zeros(1101,1);
        ptemp=zeros(1101,length(BinA));
        for i=1:1101
        El=0.001*(i-1);
        ka=kmol*El+kabg;kd=kmol*(1-El)+kdbg;
        if ka<0 %this is another insurance policy in addition to the above, prevent the situation:ka<0,kd<0 from happening
            ka=1;
        end
        if kd<0
            kd=1;
        end
        pa=-ka+BinA.*log(ka)-log(2*pi.*BinA)./2-BinA.*log(BinA)+BinA;
        pd=-kd+BinD.*log(kd)-log(2*pi.*BinD)./2-BinD.*log(BinD)+BinD;
        ptol=pa+pd;
        pe=exp(ptol);
        ps(:,i)=pe;
        p(i)=sum(pe); 
        ptemp(i,:)=pe;
%         Peverypoint(i,FretFramePosition(j,1):FretFramePosition(j,2))=pe; %record the poisson distribution of every data point
        end
        Peverypoint(:,FretFramePosition(j,1):FretFramePosition(j,2))=ptemp;
        [~,ps]=max(ps,[],2);
        ps=ps.*0.001;
        dif(j)=abs(mean(data2)-mean(ps)); %calculate the difference between original result and recalculated result
        if dif(j)>0.01 %if the difference is significant
            ps=data2;
            BinD=kmol*(1-ps)+kdbg; %construct a virtual intensity-frame trajectory using the original result
            BinA=kmol*ps+kabg;
            ps=zeros(length(BinA),1101);
            p=zeros(1101,1);
            if kmol*-0.1+kdbg<1
            Add=-(kmol*-0.1+kdbg)+1;
            kdbg=kdbg+Add;
            BinD=BinD+Add;
            end
            if ~isempty(find(BinA<=0, 1)) %BinA<=0 is not allowed, we add the background level to prevent this from happening
                Add=(-min(BinA))+1;
                kabg=kabg+Add;
                BinA=BinA+Add;
            end
            if ~isempty(find(BinD<=0, 1)) %BinD<=0 is not allowed, we add the background level to prevent this from happening
                Add=(-min(BinD))+1;
                kdbg=kdbg+Add;
                BinD=BinD+Add;
            end
            for i=1:1101 %recalculate FRET efficiency from this trajectory
            E=0.001*(i-1);
            ka=kmol*E+kabg;kd=kmol*(1-E)+kdbg;
            if ka<0
                ka=1;
            end
            if kd<0
                kd=1;
            end
            pa=-ka+BinA.*log(ka)-log(2*pi.*BinA)./2-BinA.*log(BinA)+BinA;
            pd=-kd+BinD.*log(kd)-log(2*pi.*BinD)./2-BinD.*log(BinD)+BinD;
            ptol=pa+pd;
            pe=exp(ptol);
            ps(:,i)=pe;
            p(i)=sum(pe);
            Peverypoint(i,FretFramePosition(j,1):FretFramePosition(j,2))=pe; %record the poisson distribution of every data point
            end
            [~,ps]=max(ps,[],2);
            ps=ps.*0.001;
            dif(j)=abs(mean(data2)-mean(ps)); %calculate the difference between original result and recalculated result again
%             dif(j)=mean(abs(data2(data3(1):data3(2),1)-ps));
%             p=p./sum(p).*(regionF(end)-regionF(1)+1);
        end
        Peverypoint(:,FretFramePosition(j,1):FretFramePosition(j,2))=Peverypoint(:,FretFramePosition(j,1):FretFramePosition(j,2))./sum(p).*length(BinA); %re-scale the distribution
        p=p./sum(p).*length(BinA);
        pall=pall+p; %record the total poisson distribution of whole trace
        for k=1:length(ps)
           yp(end+1)=ps(k); %record the FRET efficiency of the total trace
        end
        BinAA(j,1:FretFrame(j))=BinA; %record the intensity trace of Accepter channel in FRET region of every trajectory
        BinDD(j,1:FretFrame(j))=BinD; %record the intensity trace of Donor channel in FRET region of every trajectory
        Kmol(j)=kmol; %record the mean total intensity of FRET frame in both channel after correction
        Kabg(j)=kabg; %record the mean background intensity of Accepter channel
        Kdbg(j)=kdbg; %record the mean background intensity of Donor channel
    end
    
    % Ask you for another StateNumber input if StateNumber=0
    if StateNumber==0
        figure(1)
        plot(0:0.001:1.1,pall.*1000./sum(pall));
        xlim([0,1.1]);
        xlabel('Efficiency')
        ylabel('Distribution');
%         hold on;
        NeedInput=1;
        while NeedInput==1
            StateNumber=input('Input state number','s');
            IsNumber=str2double(StateNumber);
            if ~isnan(IsNumber) %StateNumber you input should be a integer
                StateNumber=str2double(StateNumber);
                if rem(StateNumber,1)==0
                    NeedInput=0;
                end
            end
        end
    end
else
    StateNumber=0; %if the folder has no .txt file, the output StateNumber will be set to 0
    return
end
end