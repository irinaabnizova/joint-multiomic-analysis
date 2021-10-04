%deal intersect DMR and DARs SMR and SARs one chromosome

clear;
display('%===========================params');
thr_vic=200000


  %==============================MAIN INPUT
  %-------------------items 
  display('items around targets (DEGs SEGs) are DARs and SARs , chr1'); 
  targetDEG= load('test_data\DEG_chr1.txt');%[dmr_ect;dmr_end;dmr_mes];%[dar_ect_uniq_1K;dar_end_uniq_1K;dar_mes_uniq_1K];%Kdmr_end;
  %--unmeth
  targetSEG= load('test_data\SEG_chr1.txt');%smrl;%sar;%dmr_mes;

  %DMR= load('test_data\DAR_chr1.txt');%[dmr_ect;dmr_end;dmr_mes];%[dar_ect_uniq_1K;dar_end_uniq_1K;dar_mes_uniq_1K];%Kdmr_end;
  %--unmeth
  %SMRL= load('test_data\SAR_chr1.txt');%smrl;%sar;%dmr_mes;

  num_item1=max(size(targetDEG));
  num_item2=max(size(targetSEG));
  total_num_deg_seg=[num_item1,num_item2]


  
 display('items are DARs lineage-specific and SARs');
 %===================target 

 DAR=load('test_data\DMR_chr1.txt');%[dar_ect_uniq_03K;dar_end_uniq_03K;dar_mes_uniq_03K];
 SAR=load('test_data\SMR_chr1.txt');%sar;%targetDAR;%pr_cg;%K27ac_Mes;

 num_tarE=max(size(DAR));
 num_tarM=max(size(SAR));
 total_num_DAR_SAR=[num_tarE,num_tarM]



%---------------------------aligned by target starts

chrN1=1;
[percE1,percM1,dis2startE1,dis2startM1,EnhTargetE1,EnhTargetM1,TargetEnhE1,TargetEnhM1,net12_1,pet12_1,ind_EnhE1,ind_EnhM1]=enhs_targets_stats_chr_sorted(chrN1,thr_vic,DAR,SAR,targetDEG,targetSEG);

display('%------------------number enh1 & enh2 around targets interacting')

dar_deg=sum(net12_1(1,:))
sar_seg=sum(net12_1(2,:))

forJacc_dmr=[num_tarE, num_item1,dar_deg(1)];
forJacc_smr=[num_tarM, num_item2,sar_seg(1)];

%disp(['Jaccard of DMR, unmeth SMRs around H3K27ac interacting, vic=',num2str(thr_vic)],'FontSize',14);
disp('Jaccard of DMR, unmeth SMRs around DARs interacting, vic=');
thr_vic
Jacc_dar_DEG=forJacc_dmr(3)/(forJacc_dmr(1)+forJacc_dmr(2)-forJacc_dmr(3))
Jacc_sar_SEG=forJacc_smr(3)/(forJacc_smr(1)+forJacc_smr(2)-forJacc_smr(3))

%==========================final stats

 
 display('percent of SARs around SEGs');
 petSAR_access=round((pet12_1(2,2)))
 display('percent of DARs around DEGs');
 petDAR=round((pet12_1(1,2)))

 
 %----------------------------------plots
 %=========================histograms in each 50bp
 figure('units','normalized','outerposition',[0 0 1 1]); 
 subplot(1,2,2);
 [bin5E,hist5E]=smooth_histo_TSS_fig(thr_vic,dis2startE1);
 title(['DARs around DEGs, vic=',num2str(thr_vic)],'FontSize',14);
 ylim([0 10]);
 subplot(1,2,1);
 [bin5M,hist5M]=smooth_histo_TSS_fig(thr_vic,dis2startM1);
 title(['accessible SARs around SEGs,vic=',num2str(thr_vic)],'FontSize',14);
 ylim([0 10]);

 %___________histograms in each 5K bp
 num_geneS=length(targetSEG(:,1));
 ngwS=num_geneS;
 [binS, count_binS,NcbS,NcbvS]=histo5K_normVicNov(thr_vic,dis2startM1,num_geneS,ngwS);

 num_geneD=length(targetDEG(:,1));
 ngwD=num_geneD;
 [binD, count_binD,NcbD,NcbvD]=histo5K_normVicNov(thr_vic,dis2startE1,num_geneD,ngwD);

 [ bin3S,hist3S ] = smooth_win3_float( count_binS,binS );
 [ bin3D,hist3D ] = smooth_win3_float( count_binD,binD );
  
 figure('units','normalized','outerposition',[0 0 1 1]); 
 subplot(1,2,2);
 plot(bin3D,hist3D,'b');%title('Normalised frequency of SARs arouns SEGs');
 title(['DARs around DEGs, vic=',num2str(thr_vic)],'FontSize',14);
 %ylim([0 1]);
 subplot(1,2,1);
 plot(bin3S,hist3S,'r');%title('Normalised frequency of SARs arouns SEGs');
 title(['accessible SARs around SEGs,vic=',num2str(thr_vic)],'FontSize',14);
 %ylim([0 10]);

