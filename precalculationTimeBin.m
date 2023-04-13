function [betaA,betaD,IbetaA,IbetaD]=precalculationTimeBin(regionF,regionC,regionB,s)
% This function is a cell function, and is a variation of the
% [...]=precalculation(...), the only difference between them is that the
% previous function is worked on photon arrival time trace, and should
% combine used with MIDA, but this function is worked on equal time binning
% trace, expecially the tirf data.Use this function to calculate the 4
% parameters that will be used in [...]=postcalculationTimeBin(...). Note
% that the input data s(:,2) contains the accepter channel data, and s(:,1)
% contains the donor channel data. regionF is the FRET frame region in both
% channels, regionC is the CrossTalk frame region in both channels. regionB
% is the background frame region in both channels. This function should be
% combined with postcalculationTimeBin(...) and Ecalculation2ndEdition.

% Manufactured by ChenTing,2020.9.16. Last update:2020.9.17

s1=s(:,1); %read two channel data, s1 is donor channel data and s2 is accepter channel data
s2=s(:,2);
bkga=sum(s2(regionB(1):regionB(end)))/(regionB(end)-regionB(1)+1); %calculate the background of each channel
bkgd=sum(s1(regionB(1):regionB(end)))/(regionB(end)-regionB(1)+1);
Ia=sum(s2(regionF(1):regionF(end)))/(regionF(end)-regionF(1)+1); %calculate the intensity in FRET region of each channel
Id=sum(s1(regionF(1):regionF(end)))/(regionF(end)-regionF(1)+1);
IbetaD=sum(s1(regionC(1):regionC(end)))/(regionC(end)-regionC(1)+1); %calculate IbetaD
Ica=sum(s2(regionC(1):regionC(end)))/(regionC(end)-regionC(1)+1);%calculate the intensity in CrossTalk region of accepter channel
Xd=(Ica-bkga)/(IbetaD-bkgd); %calculate the parameter of CrossTalk
if Xd<0 %if CrossTalk parameter is less than 0, alert user and set it to 0
    Across=(Ica-bkga)/bkga;
    Across=num2str(Across);
    disp(['The cross talk Xd is less than zero, Reference number Across is ',Across]);
    Xd=0;
end
P=(IbetaD-Id)/(Ia-Ica); %calculate intermediate parameter P
Ia0=((IbetaD-bkgd)+Xd*(IbetaD-bkgd)*P)/P; %calculate the acctpter real intensity when E=1
IbetaA=Ia0+bkga; %calculate IbetaA
betaD=IbetaD/bkgd; %calcluate betaD
betaA=IbetaA/(bkga+Xd*(IbetaD-bkgd)); %calcluate betaA
end
