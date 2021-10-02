%deal intersect DMR and DARs SMR and SARs one chromosome

clear;
display('%===========================params');
thr_vic=10000


  %==============================MAIN INPUT
  %-------------------items 
  display('items around targets (DMRs SMRs unme) are DARs and SARs , chr1'); 
  DMR= load('test_data\DAR_chr1.txt');%[dmr_ect;dmr_end;dmr_mes];%[dar_ect_uniq_1K;dar_end_uniq_1K;dar_mes_uniq_1K];%Kdmr_end;
  %--unmeth
  SMRL= load('test_data\SAR_chr1.txt');%smrl;%sar;%dmr_mes;

  num_item1=max(size(DMR));
  num_item2=max(size(SMRL));
  total_num_dmr_smr=[num_item1,num_item2]


  
 display('targets are DARs lineage-specific and SARs');
 %===================target 

 targetDAR=load('test_data\DMR_chr1.txt');%[dar_ect_uniq_03K;dar_end_uniq_03K;dar_mes_uniq_03K];
 targetSAR=load('test_data\SMR_chr1.txt');%sar;%targetDAR;%pr_cg;%K27ac_Mes;

 num_tarE=max(size(targetDAR));
 num_tarM=max(size(targetSAR));
 total_num_tarDAR_tarSAR=[num_tarE,num_tarM]



%---------------------------aligned by target starts

chrN1=1;
[percE1,percM1,dis2startE1,dis2startM1,EnhTargetE1,EnhTargetM1,TargetEnhE1,TargetEnhM1,net12_1,pet12_1,ind_EnhE1,ind_EnhM1]=enhs_targets_stats_chr_sorted(chrN1,thr_vic,DMR,SMRL,targetDAR,targetSAR);

display('%------------------number enh1 & enh2 around targets interacting')

dmr_enh=sum(net12_1(1,:))
smr_enh=sum(net12_1(2,:))

forJacc_dmr=[num_tarE, num_item1,dmr_enh(1)];
forJacc_smr=[num_tarM, num_item2,smr_enh(1)];

%disp(['Jaccard of DMR, unmeth SMRs around H3K27ac interacting, vic=',num2str(thr_vic)],'FontSize',14);
disp('Jaccard of DMR, unmeth SMRs around DARs interacting, vic=');
thr_vic
Jacc_dmr_DARs=forJacc_dmr(3)/(forJacc_dmr(1)+forJacc_dmr(2)-forJacc_dmr(3))
Jacc_smr_SARs=forJacc_smr(3)/(forJacc_smr(1)+forJacc_smr(2)-forJacc_smr(3))

%==========================final stats

 
 display('percent of unmeth SMRs around SARs');
 petSMR_unmeth=round((pet12_1(2,2)))
 display('percent of DMRs around DARs');
 petDMR=round((pet12_1(1,2)))

 
 %----------------------------------plots
 %=========================histos
 figure('units','normalized','outerposition',[0 0 1 1]); 
 subplot(1,2,2);
 [bin5E,hist5E]=smooth_histo_TSS_fig(thr_vic,dis2startE1);
 title(['DMRs around DARs, vic=',num2str(thr_vic)],'FontSize',14);
 %ylim([0 10]);
 subplot(1,2,1);
 [bin5M,hist5M]=smooth_histo_TSS_fig(thr_vic,dis2startM1);
 title(['unmethylated SMRs around SARs,vic=',num2str(thr_vic)],'FontSize',14);
 %ylim([0 10]);

%===============save

%fpm=load('C:\Users\ia1\Documents\R_work\older_figs\fig3_Jan\DMR_mes_H3k27acMes_20K_data.txt');
%C:\Users\ia1\Documents\R_work\fig3_Jan\dare_H3k27acEct_20K_data.txt
%fpm=fopen('SARs\SARs_prom_nCG_5K_data.txt','w');
%for i=1:length(bin5M),
%fprintf(fpm,'%d\t%.2f\n',bin5M(i),hist5M(i));
%end
%fclose(fpm);

%fpe=fopen('DARs\DARs_prom_nCG_5K_data.txt','w');
%for i=1:length(bin5E),
%fprintf(fpe,'%d\t%.2f\n',bin5E(i),hist5E(i));
%end
%fclose(fpe);