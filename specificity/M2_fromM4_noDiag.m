function [DM,me_med_sd_DEGDEG,me_med_sd_DEGSEG]=M2_fromM4_noDiag(Mpp)
%---------------------------make 2x2 matrix from 4x4, 15 June
%prom-enh
vec_DEG_DEG=[Mpp(1,2),Mpp(1,3),Mpp(2,1),Mpp(2,3),Mpp(3,1),Mpp(3,2)];
%prom-prom- no diagonals
%vec_DEG=[Mpp(1,2),Mpp(1,3),Mpp(2,1),Mpp(2,3),Mpp(3,1),Mpp(3,2)]
DM(1,1)=mean(vec_DEG_DEG);
DM(2,2)=Mpp(4,4);
me_med_sd_DEGDEG=[mean(vec_DEG_DEG),median(vec_DEG_DEG),std(vec_DEG_DEG)];

vec_DEG_SEG=[Mpp(1,4),Mpp(2,4),Mpp(3,4),Mpp(4,1),Mpp(4,2),Mpp(4,3)];
DM(1,2)=mean(vec_DEG_SEG);
me_med_sd_DEGSEG=[mean(vec_DEG_SEG),median(vec_DEG_SEG),std(vec_DEG_SEG)];

DM(2,1)=DM(1,2);
%DM
