function [exx,rank_perm]=put_Zero_rank_uniform_bw_a_b(a,b,N,rank_permuted)
%---------------witout replacement

%-uniform in the interval between int a and b (a,b), N samples
%m=277;
exx=ones(b,1);
%a=1;
%b=30;
%N=100;
%N=numZ;
%N=2
%a=1;
%b=m;
% get N positions and make them Zero
rank_perm=rank_permuted;
for i=1:N,
x(i)=round(a+(b-a)*rand);
exx(x(i))=0;
rank_perm(x(i))=0;
end
s=sum(exx);
%figure;
%hist(x,100)
%title('unif distr at interval');

