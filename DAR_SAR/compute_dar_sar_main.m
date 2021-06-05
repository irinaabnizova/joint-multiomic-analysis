function [DAR,SARH]=compute_dar_sar_main(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR,thr_cov,thr_dif,thr_sim,thr_accH,thr_accL)
%compute DAR SAR from counts, one chrN

%-Input= chrN and chr_name, both match 13 for this test data
%Input default parameters/thresholds: thr_cov,thr_dif,thr_sim,thr_accH,thr_accL

DAR=[];
SARH=[];% sars we are interested in are Highly accssible

%--------------------------parameters
if nargin < 13,
   display('DEfault Parameters of DAR SAR computing from GC counts in a 100 bp window, with separate access High and Low thr');
   thr_cov=25
   thr_dif=0.51;
   thr_sim=0.36; % for pu< than the_sim
   thr_accH=0.35;
   thr_accL=0.20;
end


display('Variables and parameters');

%===========================variables to define manually here
   chrN
   name_chr
   thr_cov
   thr_accHL=[thr_accH,thr_accL]
   thr_sim_dif=[thr_sim,thr_dif]
%--------------------folders and names for input and output
   
   %winFolder=('test_data\windows\');
   %my_folder='test_data\counts\';
   %folder='test_data\output_DARSAR\';
    
   %FilenameWin=sprintf('wine100_%s.txt',name_chr);
   %FilenameEct  = sprintf('me_unme_nposm_nposu_Ect1_%s.txt',name_chr);
   %FilenameEnd  = sprintf('me_unme_nposm_nposu_End5_%s.txt',name_chr);
   %FilenameMes  = sprintf('me_unme_nposm_nposu_Mes_%s.txt',name_chr);
   
   %outDAR=sprintf('DAR_EctEndMes_c25_Pu051_%s.txt',name_chr);
   %outSAR=sprintf('SARH_EctEndMes_c25_Pu036_%s.txt',name_chr);
    
%--------------------params

display('1. Filter precomputed counts by a given coverage');
       
       %file_win=fullfile(winFolder, FilenameWin);
       %file_countEct=fullfile(my_folder, FilenameEct);
       %f%ile_countEnd=fullfile(my_folder, FilenameEnd);
       %file_countMes=fullfile(my_folder, FilenameMes);
       
    ew=fopen(file_win);
    if ew<0,
        display('wrong/non existent input1 window file');
        return
    end
    e_ect=fopen(file_countEct);
    if e_ect<0,
        display('wrong/non existent input2 ect count GC file');
        return
    end
    e_end=fopen(file_countEnd);
    if e_end<0,
        display('wrong/non existent input3 end count GC file');
        return
    end
    e_mes=fopen(file_countMes);
    if e_mes<0,
        display('wrong/non existent input4 mes count GC file');
        return
    end
    
       % win=load(fullfile(winFolder, FilenameWin));
       %Y_muEc = load(fullfile(my_folder, FilenameEct));
       %file_win=fullfile(winFolder, FilenameWin);
       %file_countEct=fullfile(my_folder, FilenameEct);

      [winn,level_BEFneg,covEFneg,BE_npos_meNeg,BE_npos_unmeNeg,fracZ_lowE]=prefilter_count_chr(thr_cov,file_win,file_countEct);
      [winn,level_BEnFneg,covEnFneg,BEn_npos_meNeg,BEn_npos_unmeNeg,fracZ_lowEn]=prefilter_count_chr(thr_cov,file_win,file_countEnd);
      [winn,level_BCFneg,covCFneg,BC_npos_meNeg,BC_npos_unmeNeg,fracZ_lowM]=prefilter_count_chr(thr_cov,file_win,file_countMes);

      fracZ_low_EcEnM=[fracZ_lowE;fracZ_lowEn;fracZ_lowM]
   
display('2. DEfine DARs and  equal SARs HL with standard params') 
   %[tot, DAR, DAR_raw,SAR_raw, SARH,pu,ch,amp,SARL]=dar_SARHL_tot_lev_levGC_cov_PuHL(thr_sim,thr_dif,thr_accH,thr_accL,winn,level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg, chrN);
   [tot, DAR, DAR_raw, SAR_raw,SARH,pu,ch,amp,SARL]=dar_sarHL_levGC_dominance(thr_sim,thr_dif,thr_accH,thr_accL,winn,level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg, chrN);
 
   
   
   %------------------------OUTPUT
    %SAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %     1   2  3  4     5        6        7       8

   si_SARH=size(SARH);
   si_SARL=size(SARL);
   si_DAR=length(DAR(:,1))
 
   %===============write DARs SARs into output files, with accessibility level per lineage
    [a]=save_chromatin_layers(DAR,outDAR);
    [b]=save_chromatin_layers(SARH,outSAR);
     