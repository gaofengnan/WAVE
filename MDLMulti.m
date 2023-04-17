function [sfit,MDLN]=MDLMulti(filename) 
%This function is a variant of MDL function. Used to get states of E traces
%from a data set, using minimum description length (MDL) principle.

%Input parameter:filename is the raw data structure, every filename(i).ss
%record one data trace. Output parameter:sfit is the fit trace structure,
%sfit(i).ss(:,j) record the fit trace corresponding to raw data S(i).ss(or
%filename(i).ss) and reference number MDLN(i,j). MDLN is the reference
%number of all data trace in each merger.

%Manufactured by ChenTing,2022.11.07. Last update:2022.11.07
%2022.11.07 update: The results of change_point_detection is the frame
%before the intensity change began, the fitting curve part of the previous
%code ignored this problem. For example, frame 1-6 have intensity I1, frame
%7-9 have intensity I2. The change point is frame number 6. The previous
%fitting curve part fit frame 1-6 with I1, and then fit 6-9 with I2, which
%means the intensity of frame number 6 will be fitted with I2 arther than
%I1. This is incorrect. We fix this bug in this update, change fitting
%curve part in line number 58-69. Note that, in most situation the first
%frame will not going to be a change point, but if that happen, we should
%assign a separate value to the first frame, and assign other values to
%subsequent frames starting with the second frame. Change point frame
%number 1 is represented as -1 in transition and is treated as frame number
%2 when fitting curve. See line number 60-62.

%% load file and get the important parameters
S=filename;
L=length(S);
sd=zeros(1,L); %record noise level of every data trace
   N2=zeros(1,L);
 N=0; %N records the total number of data points after all data traces are connected
V=zeros(1,L); %record domain size of every data trace
Points=struct([]); %record change points of every data trace
for i=1:L
    s=S(i).ss;
    d=diff(s); %calculate the noise level use Haar
    d=abs(d);
    d=sort(d);
    sd(i)=d(round(0.682*length(d)))/sqrt(2); %use the white noise distribution and cumulative distribution function to get noise level
    if sd(i)==0 %used for test
        sd(i)=1;
    end
     N2(i) = length(s);
    N=N+length(s); %update N
    V(i) = max(s)-min(s); %get the domain size of this trace
    Points(i).points=change_point_detection(s); %get the change points of this trace
end
G=StategroupingMulti(S,Points); %get the state grouping G of this data set

%% get the fitted data
MDLN=zeros(L,length(G)); %record reference number of all data trace in each merger.
sfit=struct([]);
for i=1:L %initialize sfit
    sfit(i).ss=zeros(length(S(i).ss),length(G));
end

for i=1:length(G) %according to the G structure, G(i) means the i-1th merger 
    for j=1:length(G(i).g) %set the fitted data
        for k=1:(length(G(i).g(j).gg)-1)/3 %the structure of G(i).g(j).gg is [trace number1,start1,end1,trace number2,start2,end2,...,mean intensity], length:3n+1 
            if G(i).g(j).gg(3*k-1)==-1 %if frame number 1 is regarded as a change point,treated it as frame number 2 when calculating mean intensity and state length
                sfit(G(i).g(j).gg(3*k-2)).ss(2:G(i).g(j).gg(3*k),i)=G(i).g(j).gg(end);
            else %calculate first section's information, if the end frame is regarded as a change point,treated it as lastframe-1 when calculating mean intensity and state length
                if G(i).g(j).gg(3*k-1)==1||G(i).g(j).gg(3*k-1)+1>G(i).g(j).gg(3*k)
                    sfit(G(i).g(j).gg(3*k-2)).ss(G(i).g(j).gg(3*k-1):G(i).g(j).gg(3*k),i)=G(i).g(j).gg(end);
                else %normal situation
                    sfit(G(i).g(j).gg(3*k-2)).ss(G(i).g(j).gg(3*k-1)+1:G(i).g(j).gg(3*k),i)=G(i).g(j).gg(end);
                end
            end
        end
    end
    for z=1:L %calculate reference number of all data trace in this merger
        s=S(z).ss;
        Sfit=sfit(z).ss(:,i);
        Sd=sd(z);
        [lnDET, ST, ~] = G4th(Sfit, Sd); %get the 4th part of the GL
        F =1/2/Sd*sum(abs(s-Sfit)); %get the first part of MDL(i),the goodness of the fit
        GL = length(G(i).g)/2*log(1/2/pi)+length(G(i).g)*log(V(z)/Sd)+ST/2*log(N2(z))+0.5*lnDET; %get the cost of the model, use N(z)
%         GL = length(Sfit)/N*length(G(i).g)/2*log(1/2/pi)+length(Sfit)/N*length(G(i).g)*log(V(z)/Sd)+ST/2*log(N)+0.5*lnDET; %get the cost of the model, use total number of data points N
        MDLN(z,i)=F+GL; 
    end
    
end

end

%% Calculate the 4th part of the GL
function [lnDET, ST, SN] = G4th(sfit, sd)
[Z,~,IZ]=unique(sfit); %Z contains the unique number in sfit from small to large, IZ can be obtained by sfit=Z(IZ). For example, sfit=[1,1,3,3,2,2],Z=[1,2,3],IZ=[1,1,3,3,2,2]
SN=length(Z); %the state number
Tran=find(diff(IZ)~=0); %the state translation location
ST=length(Tran); %the number of state translation
SL=zeros(1,SN); 
CP=zeros(1,ST);
for i=1:SN
   SL(i)=length(IZ(IZ==i)); %the state length
end
for j=1:ST
    CP2=sfit(Tran(j))-sfit(Tran(j)+1); %the difference before and after state translation
    CP(j)=CP2^2;
end
lnDET=sum(log([SL,CP/sd^2])); %calculate the 4th part of the GL
end
