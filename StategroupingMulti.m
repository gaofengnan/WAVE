function G=StategroupingMulti(S,Points)
%This function is used to group the states for data set rather than a
%single data trace, using log likelihood merit. Called in MDLMulti.

%S is the raw data structure, every S(i).ss record one data trace, and
%Points is the change points structure, every Points(i).points record a set
%of change points found in one data trace recorded by S(i).ss.

%The output parameter G is a structure, G(i) means the the i-1th merger,
%contains TotalChangePoints-1 states, in other word, length(G(i).g) =
%TotalChangePoints-1. TotalChangePoints is the number of change points from
%all data traces.Every G(i).g(j) represent the j state under i-1th merger,
%and G(i).g(j).gg is in this form [trace number1,start1,end1,trace
%number2,start2,end2,...,mean intensity].

%Manufactured by ChenTing,2022.11.07. Last update:2022.11.07
%2022.11.07 update: The results of change_point_detection is the frame
%before the intensity change began, the calculation of mean intensity and
%frame length in the previous code ignored this problem. For example, frame
%1-6 have intensity I1, frame 7-9 have intensity I2. The change point is
%frame number 6. The previous calculation regard frame 7-9 as 4 frame,
%since 9-6+1=4, and the mean intensity is (3*I2+I1)/4, arther than I2. This
%is incorrect. We fix this bug in this update, change mean intensity and
%frame length calculation in line number 48-67. Note that, in most
%situation the first frame will not going to be a change point, but if that
%happen, we should assign a separate value to the first frame, and assign
%other values to subsequent frames starting with the second frame. Change
%point frame number 1 is represented as -1 in transition and is treated as
%frame number 2 when calculating mean intensity and frame length. See line
%number 53-55. 

G = struct([]);
G(1).g(1).gg=1; %we do not know the total change points number, so we use G(1).g(end+1).gg=... to give initial values. This can only work when the structure has an start value
M=[]; %record the length of every states in each data trace
I=[]; %record the mean intensity of every states in each data trace
for z=1:length(S)
    s=S(z).ss; %get one data trace and its corresponding change point
    points=Points(z).points;
    n=length(s);
    group1=[1,points];
    group2=[points,n];
    if ~isempty(points)
    if points(1)==1 %if the first change point is frame number 1, set it to -1
        group1(2)=-1;
    end
    end
    len=length(group1); %the length of points means it has length(points)+1 states
    m=zeros(1,len); %the length of every states in this data trace
    Ismall=zeros(1,len);
    Ismall(1)=mean(s(1:group2(1))); 
    m(1)=group2(1);
    for i=2:len
        if group1(i)==-1 %if frame number 1 is regarded as a change point,treated it as frame number 2 when calculating mean intensity and state length
            Ismall(i)=mean(s(2:group2(i)));
            m(i)=group2(i)-2+1;
        else
            if group1(i)+1>group2(i) %if the end frame is regarded as a change point,treated it as lastframe-1 when calculating mean intensity and state length
                Ismall(i)=mean(s(group1(i):group2(i)));
                m(i)=group2(i)-group1(i)+1;
            else %normal situation
                Ismall(i)=mean(s(group1(i)+1:group2(i))); %the mean intensity of every states in this data trace, corresponding to m(i) 
                m(i)=group2(i)-group1(i);
            end
        end
    end
    I=[I,Ismall]; %update I and M
    M=[M,m];
    for i=1:len
        G(1).g(end+1).gg=[z,group1(i),group2(i),Ismall(i)]; %set the initial values of G(1)
    end
end
G(1).g(1)=[]; %delete the start value of G(1)
len=length(G(1).g); %get the number of change points from all data traces

for j=2:len
    ll=ones(length(I),length(I))*(-inf); %reset the log likelihood merit matrix
    IL=zeros(length(I),length(I)); %reset the mean intensity matrix
    for p=1:(length(I)-1)
        for q=(p+1):length(I)
            [ll(p,q),IL(p,q)]=llmerit(I(p), I(q), M(p), M(q)); %set these two matrixs, because p<q, they are all upper triangular matrix
        end
    end
    [~,b]=max(ll(:)); %find the location of the max log likelihood merit in ll, note that b is a number rather than a vector, and start from column,b=(row number)*(col-1)+row
    b=b(1);
    col = ceil(b/length(I));% get the column, ceil means round up to an integer
    row = b-(col-1)*length(I);% get the row
    G(j)=G(j-1); %set the initial values of G(j),from G(j-1)
    G(j).g(row).gg(end)=[];  %row<col, we combine col to row, and replace row, empty out col, first we delate the mean intensity of G(j).g(row).gg, the row state of j-1 merger
    for x=1:(length(G(j).g(col).gg)-1)
        G(j).g(row).gg(end+1)=G(j).g(col).gg(x); %add the G(j).g(col).gg's trace number, start and end to the end of G(j).g(row).gg
    end
    I(row)=IL(row,col); %update the I(row) mean intensity
    G(j).g(row).gg(end+1)=I(row); %add the new mean intensity to the end of G(j).g(row).gg
    G(j).g(col)=[]; %empty out G(j).g(col)
    I(col)=[]; %empty out I(col)
    M(row)=M(row)+M(col); %update the M(row) length of state
    M(col)=[]; %empty out M(col)
end
end

%% This function is used to get the log likelihood merit  for a given I1 I2 m1 m2
function [ll,I] = llmerit(I1, I2, m1, m2)
I = (I1*m1 + I2*m2)/(m1+m2); %give the mean intensity
ll = (m1+m2)*I^2-m1*I1^2-m2*I2^2; %give the log likelihood merit
end