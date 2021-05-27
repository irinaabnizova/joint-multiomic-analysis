%compute DAR SAR from counts, one chrN

display('DAR SAR Ect End Mes unique, from Counts in a window: PARAMs');
chrN=13
name_chr='chr13'
thr_cov=25

display('Standard Parameters of DAR SAR defining, with separate access High and Low thr');
thr_dif=0.51;
thr_sim=0.36; % for pu< than the_sim
thr_accH=0.35;
thr_accL=0.20;

thr_accHL=[thr_accH,thr_accL]
thr_sim_dif=[thr_sim,thr_dif]

display('1. Filter precomputed counts by standard coverage=25x');

   FilenameWin=sprintf('wine100_%s.txt',name_chr);
   winFolder=('test_data\windows\');
   my_folder='test_data\counts\';
 
   FilenameEct  = sprintf('me_unme_nposm_nposu_Ect_%s.txt',name_chr);
   [winn,level_BEFneg,covEFneg,BE_npos_meNeg,BE_npos_unmeNeg,fracZ_lowE]=prefilter_count_chrN_ext(thr_cov,FilenameWin,winFolder,my_folder,FilenameEct);
   
   FilenameEnd  = sprintf('me_unme_nposm_nposu_End5_%s.txt',name_chr);
   [winn,level_BEnFneg,covEnFneg,BEn_npos_meNeg,BEn_npos_unmeNeg,fracZ_lowEn]=prefilter_count_chrN_ext(thr_cov,FilenameWin,winFolder,my_folder,FilenameEnd,chrN);
   
   FilenameMes  = sprintf('me_unme_nposm_nposu_Mes_%s.txt',name_chr);
   [winn,level_BCFneg,covCFneg,BC_npos_meNeg,BC_npos_unmeNeg,fracZ_lowM]=prefilter_count_chrN_ext(thr_cov,FilenameWin,winFolder,my_folder,FilenameMes,chrN);

   fracZ_low_EcEnM=[fracZ_lowE;fracZ_lowEn;fracZ_lowM]
   
display('2. DEfine DARs and  equal SARs HL with standard params') 
   [tot, DAR, DAR_car, SARH,pu,ch,amp,SARL]=dar_SARHL_tot_lev_levGC_cov_PuHL(thr_sim,thr_dif,thr_accH,thr_accL,winn,level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg, chrN);
    %------------------------OUTPUT
    %SAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %     1   2  3  4     5        6        7       8

   si_SARH=size(SARH)
   si_SARL=size(SARL)
   si_DAR=size(DAR)
 
   %===================SAVE with accessibility level per lineage
   
    folder='test_data\output_DARSAR\';
    textFilename=sprintf('DAR_EctEndMes_c25_Pu051_%s.txt',name_chr);
    [chrN]=save_DAR_EctEndMes(DAR,folder,textFilename,chrN);

    textFilename=sprintf('SARH_EctEndMes_c25_Pu036_%s.txt',name_chr);
    [chrN]=save_DAR_EctEndMes(SARH,folder,textFilename,chrN);
    
    