function [feat_existP,feat_rankP,feat_scoreP,feat_exist,feat_rank,feat_score]=quartileControl_feature_permute(va1, va2, va3,va4,f_ex,ind,list,ind1);

%-------------------------inputs

%-------------------------outputs
%====================feature-vector for original list of motifs
%feat_exist=Ex;
%feat_rank=f_list1;
%feat_score=f_list1i;- weights quartile -scored



%------------control quartile values:va1, va2, va3,va4
% apply quartile strategy for inverse ranks
%-make feature_list (FL)and its index (for names of motifs)  from feature and list
%f_all=importdata('capital_feature_name_long_2.m');%alpfabetical, 277

%f_ex=f_all;
%ind=1:length(f_ex);
%list=importdata('dmr_Mes.m');
%ind1=1:length(list);


%------------------------------OUTPUT
%feat_exist=Ex;
%feat_ind=f_list1;ran=ind within feat
%feat_rank=f_list1;
%feat_score=f_list1i;

%------------------------------------------initiate
feat_exist=[];feat_rank=[];feat_score=[];
feat_existP=[];feat_rankP=[];feat_scoreP=[];

Ex=zeros(1,length(ind));
f_list1=zeros(1,length(ind));
f_list1i=zeros(1,length(ind));
m=length(ind);

%f_list1=zeros(1,length(f_ex));
%ind_f1=zeros(1,length(list));


for i=1:m,%length(f_ex),
    for j=1:length(list),
        if strcmp(f_ex(i),list(j))>0,%f_ex(i)==list(j),
            f_list1(i)=ind1(j);% rank in the enhancer/promoter list
            ind_f1(j)=ind(i);
            Ex(i)=1;
        end
    end
end


%--------------------get quartiles n1 n2 n3
n=length(list);
n1=round(n*0.25);
n2=round(n*0.5);
n3=round(n*0.75);


%2.1--------------------------original scores

for i=1:m,%length(f_list1),
    
    %------first quartile Va1
    if f_list1(i)>=1 && f_list1(i)<= n1,
        f_list1i(i)=va1;
    end
    if f_list1(i)>n1 && f_list1(i)<= n2,
        f_list1i(i)=va2;
    end
    if f_list1(i)>n2 && f_list1(i)<= n3,
        f_list1i(i)=va3;
    end
    
     if f_list1(i)>n3 && f_list1(i)<= n,
        f_list1i(i)=va4;
     end  
end

%====================feature-vector for original list of motifs
 
feat_exist=Ex;
feat_rank=f_list1;
feat_score=f_list1i;


%2===========================PERMUTE

%-------------------prmute rank list
%---------------------1. permute all rank features, f_list1

mi=1;
ma=max(f_list1);

%--------------------generates feature list of length m (as FLC)
%-------------------of randomly permuted ranks (1,ma)-unif distr
N=m;
ranks_permuted=round(uniform_bw_a_b(mi,ma,N));

%----------------randomly insert numZ zeros in this feature rank
numZ=m-ma;
N1=numZ;
%[Ex_permuted]=put_Zero_uniform_bw_a_b(1,m,round(2.5*N1));
[Ex_perm,rank_perm]=put_Zero_rank_uniform_bw_a_b(1,m,round(2.5*N1),ranks_permuted);



%-----------------------------make ranki for permuted list
%Ex_perm,rank_perm
%2.2-------------------------permuted scores
rank_permi=rank_perm;% iinitiate to make same size
for i=1:m,%length(f_list1),
    
    %------first quartile Va1
    if rank_perm(i)>=1 && rank_perm(i)<= n1,
        rank_permi(i)=va1;
    end
    if rank_perm(i)>n1 && rank_perm(i)<= n2,
        rank_permi(i)=va2;
    end
    if rank_perm(i)>n2 && rank_perm(i)<= n3,
        rank_permi(i)=va3;
    end
    
     if rank_perm(i)>n3 && rank_perm(i)<= n,
        rank_permi(i)=va4;
     end  
end
 
feat_existP=Ex_perm;
%feat_ind=f_list1;
feat_rankP=rank_perm;
feat_scoreP=rank_permi;





