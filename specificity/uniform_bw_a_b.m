function [x]=uniform_bw_a_b(a,b,N)

%-uniform in the interval between int a and b (a,b)



%clear x;
for i=1:N,
x(i)=a+(b-a)*rand;
end

%figure;
%hist(x,100)
%title('unif distr at interval');

