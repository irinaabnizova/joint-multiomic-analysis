function [ bin3,hist3 ] = smooth_win3_float( hist,bins )
%smoother in a sliding window 3bp
 bin3=[];
 hist3=[];

%get average in three bins

%-------------First bin av over first two
bin31(1)=bins(1);
bin_width=bins(2)-bins(1);

for i=2:length(hist),
    bin31(i)=bin31(i-1)+bin_width;
end

hist3(1)=(hist(1)+hist(2))/2;
%hist3(1)=(hist(1)+hist(2)+hist(3))/3;
for i=2:length(hist)-1,
    
hist3(i)=(hist(i-1)+hist(i)+hist(i+1))/3;

end

%---------------Last bin
hist3(length(hist))=(hist(length(hist))+hist(length(hist)-1))/2;

bin3=bin31(1:length(hist3));
%hist33=round(hist3);

