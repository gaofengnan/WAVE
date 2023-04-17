%NonequilibriumCompare 23.2.1
% This script is used to evaluate our calculation results, using data sets
% produced by MakeNonequilibriumHMMPoissonData
% The workflow is like follows:
% Step 1:use MakeNonequilibriumHMMPoissonData to generate a data set
% Step 2:use TestNonequilibrium3 to analysis this data set
% Step 3:use NonequilibriumCompare to evaluate the calculation results
% The output: 
%            EvaluationResult: is a column vector, like follows: 
%      [The number of Pass; The mumber of Acceptable; The number of Failed]
%            NonequilibriumCompareParameter: record the file name, ground
%            truth in the first 3 rows, as well as our calculation results
%            in last 2 rows.
% Note: This script can only be used after TestNonequilibrium3, because
% some variables in this script are inherited from TestNonequilibrium3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mode controls which evaluation criteria you want to use, the 
% corresponding relationship is as follows: 
% mode=1 :evaluate the results of data sets produced by non-overlapping
%         7-state system [A-1],non-overlapping 4-state system [A-2] in
%         MakeNonequilibriumHMMPoissonData 
% mode=2 :evaluate the results of data sets produced by partially
%         overlapping 5-state system [B] in
%         MakeNonequilibriumHMMPoissonData 
% mode=3 :evaluate the results of data sets produced by partially
%         overlapping 4-state system [C] in
%         MakeNonequilibriumHMMPoissonData 
% mode=4 :evaluate the results of data sets produced by partially
%         overlapping 4-state system [D] in 
%         MakeNonequilibriumHMMPoissonData 
mode=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the ground truth of this data set, and combine it with our
% calculation results.
DataSetInformation=load([file_path,'\DataSetInformation.mat']); 
DataSetInformation=DataSetInformation.DataSetInformation;
l=length(txt_path_list);
[~,n]=size(DataSetInformation);
if l~=n
    disp('Wrong data set')
    return
end
Seq=zeros(1,n); %Seq is used to reorder our calculation results, make them in the same order as in GT 
for i=1:n
    for j=1:n
        Name=DataSetInformation{1,i};
        if strcmp(txt_path_list(j).name,Name)
            Seq(i)=j;
            break
        elseif strcmp(txt_path_list(j).name,Name(1:end-15))
            Seq(i)=j;
            break
        elseif strcmp(txt_path_list(j).name(1:end-12),Name(1:end))
            Seq(i)=j;
            break
        end
    end
end
if ~isempty(find(Seq==0, 1))
    disp('Wrong data set')
    return
end

NonequilibriumCompareParameter=cell(5,n); %construct NonequilibriumCompareParameter by GT and calculation results
NonequilibriumCompareParameter(1:3,1:n)=DataSetInformation(1:3,1:n);
for j=1:n
    NonequilibriumCompareParameter{4,j}=AfterLengthResult(Seq(j));
    NonequilibriumCompareParameter{5,j}=StateBefore(Seq(j));
end


yfitreal=[]; %load the real fitting trajectories from GT
for i=1:n
    A=isstrprop(txt_path_list(i).name,'digit');
    B=txt_path_list(i).name(A);
    lo=str2num(B);
    yfitreal=[yfitreal;DataSetInformation{4,lo}];
end

% evaluate the calculation results, using criteria you select
Information=cell2mat(NonequilibriumCompareParameter(2:end,1:end));

if mode==1 %evaluate the results of data sets produced by non-overlapping 7-state system [A-1],non-overlapping 4-state system [A-2] in MakeNonequilibriumHMMPoissonData
    FrameDifference=abs(Information(1,:)-Information(3,:));
    Count1=2.*(FrameDifference>2&FrameDifference<=5);
    Count2=4.*(FrameDifference>5);
    StateDifference=abs(Information(2,:)-Information(4,:));
    Count3=StateDifference>=0.05;
    Count=Count1+Count2+Count3;
    EvaluationResult=[length(find(Count<2));length(find(Count==2));length(find(Count>2))];
    figure(5000)
    bar(1,EvaluationResult,0.9,'stacked','EdgeColor','k')
    disp(['Pass:', num2str(EvaluationResult(1)), ' ;Acceptable:', num2str(EvaluationResult(2)),' ;Failed:', num2str(EvaluationResult(3))])
elseif mode==2 %evaluate the results of data sets produced by partially overlapping 5-state system [B] in MakeNonequilibriumHMMPoissonData
    FrameDifference=Information(1,:)-Information(3,:);
    L1=length(find(abs(FrameDifference)<=2));
    L2=length(find(FrameDifference<-2&-5<=FrameDifference));
    L3=length(find(FrameDifference<-5));
    L4=find(FrameDifference>2);
    for i=1:length(L4)
        if Seq(L4(i))==1
            yfitrealpart=yfitreal(1:TraceLength(1));
            yfitpart=yfit(1:TraceLength(1));
            ypart=y(1:TraceLength(1));
        else
            yfitrealpart=yfitreal(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
            yfitpart=yfit(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
            ypart=y(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
        end
        yfitrealpart2=yfitrealpart(1:TraceLength(Seq(L4(i)))-Information(3,L4(i)));
        if length(find(yfitrealpart2==0.9))<=2
            L1=L1+1;
        elseif length(find(yfitrealpart2==0.9))<=5
            L2=L2+1;
        else
            L3=L3+1;
        end
        figure(i)
        for k=1:4
        plot([1,length(ypart)],[k*0.2,k*0.2],'k')
        hold on;
        end
        plot([1,length(ypart)],[0.9,0.9],'k')
        plot(ypart,'k');hold on;plot(yfitpart,'r');plot(yfitrealpart,'b');
        plot([TraceLength(Seq(L4(i)))-Information(1,L4(i))+1,TraceLength(Seq(L4(i)))-Information(1,L4(i))+1],[0,1],'b','LineWidth',1.5)
        plot([TraceLength(Seq(L4(i)))-Information(3,L4(i))+1,TraceLength(Seq(L4(i)))-Information(3,L4(i))+1],[0,1],'r','LineWidth',1.5)
        ylim([0,1])
        hold off;
    end
    figure(5000)
    EvaluationResult=[L1;L2;L3];
    bar(1,EvaluationResult,0.9,'stacked','EdgeColor','k')
    disp(['Pass:', num2str(EvaluationResult(1)), ' ;Acceptable:', num2str(EvaluationResult(2)),' ;Failed:', num2str(EvaluationResult(3))])
elseif mode==3 %evaluate the results of data sets produced by partially overlapping 4-state system [C] in MakeNonequilibriumHMMPoissonData
    FrameDifference=Information(1,:)-Information(3,:);
    L1=length(find(abs(FrameDifference)<=2));
    L2=length(find(FrameDifference>2&FrameDifference<=5)); 
    L3=length(find(FrameDifference>5));
    L4=find(FrameDifference<-2);
%     L5=find(FrameDifference>5); 

    for i=1:length(L4)
        if Seq(L4(i))==1
            yfitrealpart=yfitreal(1:TraceLength(1));
            yfitpart=yfit(1:TraceLength(1));
            ypart=y(1:TraceLength(1));
        else
            yfitrealpart=yfitreal(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
            yfitpart=yfit(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
            ypart=y(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
        end
        yfitrealpart2=yfitrealpart(TraceLength(Seq(L4(i)))-Information(3,L4(i))+1:end);
        if length(find(yfitrealpart2==0.8|yfitrealpart2==0.4))<=2
            L1=L1+1;
        elseif length(find(yfitrealpart2==0.8|yfitrealpart2==0.4))<=5
            L2=L2+1;
        else
            L3=L3+1;
%             figure(i)
%             plot(ypart,'k');hold on;plot(yfitpart,'r');plot(yfitrealpart,'b');
%             plot([TraceLength(Seq(L4(i)))-Information(1,L4(i))+1,TraceLength(Seq(L4(i)))-Information(1,L4(i))+1],[0,1],'b','LineWidth',1.5)
%             plot([TraceLength(Seq(L4(i)))-Information(3,L4(i))+1,TraceLength(Seq(L4(i)))-Information(3,L4(i))+1],[0,1],'r','LineWidth',1.5)
%             ylim([0,1])
%             hold off;
        end

    end

    figure(5000)
    EvaluationResult=[L1;L2;L3];
    bar(1,EvaluationResult,0.9,'stacked','EdgeColor','k')
    disp(['Pass:', num2str(EvaluationResult(1)), ' ;Acceptable:', num2str(EvaluationResult(2)),' ;Failed:', num2str(EvaluationResult(3))])
elseif mode==4 %evaluate the results of data sets produced by partially overlapping 4-state system [D] in MakeNonequilibriumHMMPoissonData
    FrameDifference=Information(1,:)-Information(3,:);
    L1=length(find(abs(FrameDifference)<=2));
    L2=length(find(FrameDifference<-2&-5<=FrameDifference));
    L3=length(find(FrameDifference<-5));
    L4=find(FrameDifference>2);
    for i=1:length(L4)
        if Seq(L4(i))==1
            yfitrealpart=yfitreal(1:TraceLength(1));
            yfitpart=yfit(1:TraceLength(1));
            ypart=y(1:TraceLength(1));
        else
            yfitrealpart=yfitreal(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
            yfitpart=yfit(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
            ypart=y(sum(TraceLength(1:Seq(L4(i))-1))+1:sum(TraceLength(1:Seq(L4(i)))));
        end
        yfitrealpart2=yfitrealpart(1:TraceLength(Seq(L4(i)))-Information(3,L4(i)));
        if length(find(yfitrealpart2==0.8|yfitrealpart2==0.2))<=2
            L1=L1+1;
        elseif length(find(yfitrealpart2==0.8|yfitrealpart2==0.2))<=5
            L2=L2+1;
        else
            L3=L3+1;
        end
%         figure(i)
%         plot(ypart,'k');hold on;plot(yfitpart,'r');plot(yfitrealpart,'b');
%         plot([TraceLength(Seq(L4(i)))-Information(1,L4(i))+1,TraceLength(Seq(L4(i)))-Information(1,L4(i))+1],[0,1],'b','LineWidth',1.5)
%         plot([TraceLength(Seq(L4(i)))-Information(3,L4(i))+1,TraceLength(Seq(L4(i)))-Information(3,L4(i))+1],[0,1],'r','LineWidth',1.5)
%         ylim([0,1])
%         hold off;
    end
    figure(5000)
    EvaluationResult=[L1;L2;L3];
    bar(1,EvaluationResult,0.9,'stacked','EdgeColor','k')
    disp(['Pass:', num2str(EvaluationResult(1)), ' ;Acceptable:', num2str(EvaluationResult(2)),' ;Failed:', num2str(EvaluationResult(3))])
end