function [percE,percM,dis2startE,dis2startM,EnhTargetE,EnhTargetM,TargetEnhE,TargetEnhM,net12,pet12,ind_EnhE,ind_EnhM,npg_12]=enhs_targets_stats_chr_sorted(chrN,thr_vic,DARE4_100,DARC4_100,targetE,targetM)

 %-------------------dars -targets for one chrom (sorted now 1 oct 2021 within current chromosome)
 
 %chrN=1
 %num_enh_tar_12=[num_int_Enh_tarE;num_int_Enh_tarM]
 %perc_enh_tar_12=[per_int_Enh_tarE;per_int_Enh_tarM]

 %net12=num_enh_tar_12;
 %pet12=perc_enh_tar_12;

 %-----------------get current chromosome
 [DARC4_100_chr1,indM]=separate_chr(DARC4_100, chrN);
 [DARE4_100_chr1,indE]=separate_chr(DARE4_100, chrN);
%-------------------sort it!!!!
 DARC4_100_chr=sortrows(DARC4_100_chr1,2);
 DARE4_100_chr=sortrows(DARE4_100_chr1,2);
 
 %1 3149333 400 0.35 0.23 40 31 1.740000e+00 2375
 %1 3160133 400 0.16 0.45 31 29 -2.940000e+00 -3442
 %1 3173833 300 0.25 0.00 28 22 9.300000e-01 2583
 
  num_items_12=[length(indE),length(indM)];

  n=length(DARE4_100_chr(1,:));
  dare=DARE4_100_chr(:,2:n);
  darc=DARC4_100_chr(:,2:n);
 
  %0.2--------------------- targets per chromosome

  [target_chrE1,ind_targetE]=separate_chr(targetE,chrN);%pos_enhM,pos_enhEc,pos_enhEn);
  [target_chrM1,ind_targetM]=separate_chr(targetM,chrN);%pos_enhM,pos_enhEc,pos_enhEn);
  %--------------sort
   target_chrE=sortrows(target_chrE1,2);
   target_chrM=sortrows(target_chrM1,2);
   %         7     4652001     4652500
   %         7     6744001     6745000
   %         7     6824001     6825000
   %         7    16708001    16708500
  

%==================================main annotation function

targetEE=target_chrE;
targetMM=target_chrM;
ne=length(target_chrE(:,1));
nm=length(target_chrM(:,1));

%display('original amount of targets, type1 vs type2');
num_tar_12=[ne,nm];
ie=(1:ne)';im=(1:nm)';

%------------------------------with indexes per chromosome
targetEE=[targetEE ie];
targetMM=[targetMM im];
add=0;%diatsnce from the edges
% E for dare
[count_Enh_targetE,EnhTargetE,TargetEnhE,num_int_Enh_tarE,per_int_Enh_tarE,dis2startE,percE,ind_EnhE]=around_target_dis2start(thr_vic,dare,targetEE,add);

%M for darc
[count_Enh_targetM,EnhTargetM,TargetEnhM,num_int_Enh_tarM,per_int_Enh_tarM,dis2startM,percM,ind_EnhM]=around_target_dis2start(thr_vic,darc,targetMM,add);


%=======================

num_enh_tar_12=[num_int_Enh_tarE;num_int_Enh_tarM];

perc_enh_tar_12=[per_int_Enh_tarE;per_int_Enh_tarM];

net12=num_enh_tar_12;
pet12=perc_enh_tar_12;

%-relative number per gene
if num_enh_tar_12(1,2)>0,
npg1=num_enh_tar_12(1,1)/num_enh_tar_12(1,2);
else
 npg1=0;
end

if num_enh_tar_12(2,2)>0,
npg2=num_enh_tar_12(2,1)/num_enh_tar_12(2,2);
else
 npg2=0;
end

npg_12=[npg1,npg2];

 
