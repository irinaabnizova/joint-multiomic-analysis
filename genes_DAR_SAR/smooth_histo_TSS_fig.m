function [bin5,hist5]=smooth_histo_TSS_fig(thr_vic,dis2start)
%----no figure sign, just plot
vec=dis2start;
%======================make histo, bins 50 bp
bin=[];
count=[];
bin=-thr_vic:50:(thr_vic+50);
for j=1:length(bin)-1,
    count(j)=0;
end


%------------
for i=1:length(vec),
    for j=1:length(count),
        
        if vec(i)>= bin(j) && vec(i)<bin(j+1),
             count(j)= count(j)+1;
        end
    end
end

%figure;plot(bin(1:length(bin)-1),count);hold on;

[ bin3,hist3 ] = smooth_win3_float( count,bin(1:length(bin)-1) );
[ bin4,hist4 ] = smooth_win3_float( hist3,bin3 );
[ bin5,hist5 ] = smooth_win3_float( hist4,bin4);
plot(bin(1:length(bin)-1),count,'c','LineWidth',0.5);hold on;
plot(bin5,hist5,'k','LineWidth',2);grid;
title('items frequency around targets in which vicinity they cluster');
