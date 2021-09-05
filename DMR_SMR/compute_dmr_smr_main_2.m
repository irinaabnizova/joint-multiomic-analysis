function [DMR,DMR_ect,DMR_end,DMR_mes,SMRH,SMRL]=compute_dmr_smr_main_2(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,thr_cov,thr_dif,thr_sim,thr_meH,thr_meL)
%compute DMR_lineage SMRL from counts, one chrN
%------------removed out for a while- 4 Sept
%-Input= chrN and chr_name, both match 13 for this test data
%Input default parameters/thresholds: thr_cov,thr_dif,thr_sim,thr_accH,thr_accL


%--------OUTPUTs
DMR=[];
SMRL=[];% SMRs we are interested in are Low methed

%--------------------------parameters
if nargin < 13,
   display('Default Parameters of DMR SMR computing from CG counts in a 300 bp window, with separate access High and Low thr');
   thr_cov=3
   thr_dif=0.51;
   thr_sim=0.36; % for pu< than the_sim
   thr_meH=0.7;
   thr_meL=0.5;
end


display('Variables and parameters');

%===========================variables to define manually here
   chrN
   name_chr
   thr_cov
   thr_meHL=[thr_meH,thr_meL]
   thr_sim_dif=[thr_sim,thr_dif]
%--------------------folders and names for input and output
   
   %winFolder=('test_data\windows\');
   %my_folder='test_data\counts\';
   %folder='test_data\output_DMRSMR\';
    
   %FilenameWin=sprintf('wine100_%s.txt',name_chr);
   %FilenameEct  = sprintf('me_unme_nposm_nposu_Ect6_%s.txt',name_chr);
   %FilenameEnd  = sprintf('me_unme_nposm_nposu_End6_%s.txt',name_chr);
   %FilenameMes  = sprintf('me_unme_nposm_nposu_Mes6_%s.txt',name_chr);
   
   %outDMR=sprintf('DMR_EctEndMes_c25_Pu051_%s.txt',name_chr);
   %outSMR=sprintf('SMRH_EctEndMes_c25_Pu036_%s.txt',name_chr);
    
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
        display('wrong/non existent input2 ect count CG file');
        return
    end
    e_end=fopen(file_countEnd);
    if e_end<0,
        display('wrong/non existent input3 end count CG file');
        return
    end
    e_mes=fopen(file_countMes);
    if e_mes<0,
        display('wrong/non existent input4 mes count CG file');
        return
    end
    
       % win=load(fullfile(winFolder, FilenameWin));
       %Y_muEc = load(fullfile(my_folder, FilenameEct));
       %file_win=fullfile(winFolder, FilenameWin);
       %file_countEct=fullfile(my_folder, FilenameEct);

      [winn,level_BEFneg,covEFneg,BE_npos_meNeg,BE_npos_unmeNeg,fracZ_lowE]=prefilter_count_chr(thr_cov,file_win,file_countEct);
      [winn,level_BEnFneg,covEnFneg,BEn_npos_meNeg,BEn_npos_unmeNeg,fracZ_lowEn]=prefilter_count_chr((thr_cov-1),file_win,file_countEnd);
      [winn,level_BCFneg,covCFneg,BC_npos_meNeg,BC_npos_unmeNeg,fracZ_lowM]=prefilter_count_chr(thr_cov,file_win,file_countMes);

      fracZ_low_EcEnM=[fracZ_lowE;fracZ_lowEn;fracZ_lowM]
   
      display('2. DEfine DMRs and  equal SMRs HL with standard params') 

     [tot1,SMR,SMRH,SMRL,DMR,DMR_ect,DMR_end,DMR_mes,pu,ch,amp]=DMR_SMR_lev_cov_levCG_PuHL(thr_dif,thr_sim,thr_meH,thr_meL,winn,level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg);

 
   
     display('3.---------write DMRs SMRs into output files, with accessibility level per lineage');
   
     [chrN]=save_DMR_SMR_EctEndMes_H07_L05(DMR_ect,DMR_end,DMR_mes,SMRH,SMRL,name_chr,chrN)

   
     