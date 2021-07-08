function [M2pp,rel_err,DSP,errP,prob]=permute_original_loop(N,va1, va2, va3,va4,FLC,ind,listS,indS,list1,ind1,list2,ind2,list3,ind3)

DSP=[];
errP=[];

for i=1:N,
%S-------------SEG proms
[feat_existPS,feat_rankPS,feat_scorePS,feat_existS,feat_rankS,feat_scoreS]=quartileControl_feature_permute(va1, va2, va3,va4,FLC,ind,listS,indS);
%1--------------------------p Ect
[feat_existP1,feat_rankP1,feat_scoreP1,feat_exist1,feat_rank1,feat_score1]=quartileControl_feature_permute(va1, va2, va3,va4,FLC,ind,list1,ind1);
%2--------------------------p End
[feat_existP2,feat_rankP2,feat_scoreP2,feat_exist2,feat_rank2,feat_score2]=quartileControl_feature_permute(va1, va2, va3,va4,FLC,ind,list2,ind2);
%3----------------------------------p Mes
[feat_existP3,feat_rankP3,feat_scoreP3,feat_exist3,feat_rank3,feat_score3]=quartileControl_feature_permute(va1, va2, va3,va4,FLC,ind,list3,ind3);

%3--------------------sim matrix
%display('similarity bw SEGs and DEGs prom-prom Before GE thr: Ect End Mes  SEG');

%format short;
[Mpp]=sim_matrix_sameType_1(feat_scoreS,feat_score1,feat_score2,feat_score3);
[M2pp,mms_dd,mms_ds]=M2_fromM4_noDiag(Mpp);
rel_err=(M2pp(1,2)-M2pp(1,1))/M2pp(1,1);

%display('similarity bw SEGs and DEGs prom-prom Permuted: Ect End Mes  SEG');

[MppP]=sim_matrix_sameType_1(feat_scorePS,feat_scoreP1,feat_scoreP2,feat_scoreP3);
[M2ppP,mms_dd,mms_ds]=M2_fromM4_noDiag(MppP);
rel_errP=(M2ppP(1,2)-M2ppP(1,1))/M2ppP(1,1);

DSP(i,:)=[M2ppP(1,1),M2ppP(1,2)];
errP(i)=rel_errP;
%err(i)=rel_err;
end
%err=[errP;err'];
%figure;histogram(errP,30);

%---------------------compute number of errP >= rel_err
nm=0;
for i=1:N,
    if errP(i)>=rel_err,
        nm=nm+1;
    end
end
prob=nm/N;
    
