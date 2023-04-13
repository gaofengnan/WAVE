function points = change_point_detection(eff)
%This function uses student-t test to find the change points in eff, result
%in points, the locations of change points

%eff is in the data form and contains one column, and this function is the call function form.

%Note that, points is a row vector

%Manufactured by Bo Shuang Landes,Research Group Rice University,Department of Chemistry,July 2014.Modified by ChenTing,2020.8.21. Last update:2020.8.21

%% The main function is to detect all the change points in a trace
y = abs(diff(eff));
y = sort(y);
sd = y(round(0.682*numel(diff(eff))))/sqrt(2); % calculate the noise level
points = recursion1(eff,sd,[], 0);% recursively detect all the change points
points = sort(points);
end

%% calculate the changepoints
function points = recursion1(eff, sd, points, counter)
%{
tau998_table = [1:30, 40, 50, 60, 80, 100, 1000;...
    318.31, 22.327, 10.215, 7.173, 5.893, 5.208, 4.785, 4.501, 4.297, 4.144,...
    4.025, 3.93, 3.852, 3.787, 3.733, 3.686, 3.646, 3.61, 3.579, 3.552,...
    3.527, 3.505, 3.485, 3.467, 3.45, 3.435, 3.421, 3.408, 3.396, 3.385,...
    3.307, 3.261, 3.232, 3.195, 3.174, 3.098];% the t-distribution
%}
N = numel(eff);
tau998 = 3.733;
%tau998 = 1.5;
if N < 2% only one point left in the segment, stop searching for the change point
    return
else
    llr = change_point_wavelet(eff,sd);
    [Z, k] = max(abs(llr));
    if Z > tau998
        counter = 0;
        points(end+1) = k;
        points1 = recursion1(eff(1:k), sd, [], counter); %break the data into two segments at the changepoint, and find the changepoint in two segments
        points2 = recursion1(eff(k+1:end), sd, [], counter);
        points = [points, points1, points2+k];
    elseif counter < 3 %if Z>3.787 is not satisfied, we need to divide the data into two segments in the middle, and find the changepoint in two segments
        counter = counter +1; %This process will repeat 3 times if Z>3.787 is not satisfied after we divide the data into two segments
        k = floor(numel(eff)/2); %divide the data into two segments in the middle
        points1 = recursion1(eff(1:k), sd, [], counter);
        points2 = recursion1(eff(k+1:end), sd, [], counter);
        points = [points, points1, points2+k];
    else
        return
    end
end
end

%% This function uses the StudentTtest to give the differences in the given data-eff, with the noise level-sd
function llr = change_point_wavelet(eff, sd)
N = numel(eff);
llr = zeros(size(eff));
for i = 1 : N-1
    I1 = mean(eff(1:i));
    I2 = mean(eff(i+1:end));
    llr(i) = (I2 - I1)/sd/sqrt(1/i+1/(N-i));
end
end