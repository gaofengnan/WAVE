function [J,E]=postcalculationTimeBin(s,betaA,betaD,IbetaA,IbetaD,regionF)
% This function is a cell function, and is used to calculate the FRET
% efficiency and the Fisher's information of each frame in FRET region.
% Must be combined with [...]=precalculationTimeBin(...) and
% Ecalculation2ndEdition. Output J is the Fisher's information of FRET
% region frame and E is the FRET efficiency of FRET region frame, both in
% the order of frame number.

% Manufactured by ChenTing,2020.9.16. Last update:2020.9.17

s1=s(:,1); %read data from two channel, s2 is Accepter channel data and s1 is Donor channel data
s2=s(:,2); 
detaT=1; %actually, the real detaT is 0.1s, but in precalculationTimeBin, we only use mean function to calculate the intensity, so the real intensity is 10 times larger. Because of this reason, detaT is set to 1s
J=zeros(1,length(regionF));
E=zeros(1,length(regionF));
for i=regionF(1):regionF(end) %calculate the Fisher's information and FRET efficiency of every frame in FRET region
    na=s2(i);
    nd=s1(i);
    Et=(IbetaD*na-IbetaA*nd*1/betaA)/(IbetaD*na*(1-1/betaD)+IbetaA*nd*(1-1/betaA));
    J(i+1-regionF(1))=IbetaD*detaT*(1-1/betaD)^2/(1-Et*(1-1/betaD))+IbetaA*detaT*(1-1/betaA)^2/(Et*(1-1/betaA)+1/betaA); %calculate j
    E(i+1-regionF(1))=Et;
end
end