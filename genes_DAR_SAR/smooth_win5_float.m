function [ bin5,hist5 ] = smooth_win5_float( hist,bins )
%UNTITLED3 shift min insert is indirectly taken by appropriate bins input
%   Detailed explanation goes here

 bin5=[];
 hist5=[];

%get average in three bins

%-------------First bin av over first two
bin51(1)=bins(1);
bin_width=bins(2)-bins(1);

for i=2:length(hist),
    bin51(i)=bin51(i-1)+bin_width;
end

hist5(1)=(hist(1)+hist(2))/2;

%-------------------Second bin: av over three

hist5(2)=(hist(1)+hist(2)+hist(3))/3;

%-------------------Third bin: av over four

hist5(3)=(hist(1)+hist(2)+hist(3)+hist(4))/4;


%-----------av over 5
for i=4:length(hist)-4,
    
hist5(i)=(hist(i-1)+hist(i)+hist(i+1)+hist(i+2)+hist(i+3))/5;

end

%---------------Last bins
hist5(length(hist))=(hist(length(hist))+hist(length(hist)-1))/2;
hist5(length(hist)-1)=(hist(length(hist))+hist(length(hist)-1)+hist(length(hist)-2))/3;
hist5(length(hist)-2)=(hist(length(hist))+hist(length(hist)-1)+hist(length(hist)-2)+hist(length(hist)-3))/4;
hist5(length(hist)-3)=(hist(length(hist))+hist(length(hist)-1)+hist(length(hist)-2)+hist(length(hist)-3)+hist(length(hist)-4))/5;

bin5=bin51(1:length(hist5));
%hist53=round(hist5);

