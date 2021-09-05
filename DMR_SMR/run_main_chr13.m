%deal from counts Ect

%deal  get DMR SMR from counts chrN

display('DMR SMR Ect End Mes unique, Count Parameters');
chrN=13
name_chr='chr13'
thr_cov=3

display('Standard Parameters of SMR defining, with separate (1-me) High and Low thr');

%thr_dif=0.5
%thr_meL=0.35
%thr_meH=0.6
%thr_sim=0.35

thr_sim=0.35; % for pu< than the_sim; pu is computed for (1-me)
thr_dif=0.5;% for pu>thr-dif



%======================meth level : not too sharp, sadly....
thr_meH=0.7;%
thr_meL=0.5;

thr_DS_pu=[thr_dif,thr_sim]

 FilenameWin=sprintf('wine100_%s.txt',name_chr);
 winFolder=('test_data\windows\');
 file_win=(fullfile(winFolder, FilenameWin));
 
 my_folder='test_data\counts\';

FilenameEct  = sprintf('me_unme_nposm_nposu_meth_Ect6_%s.txt',name_chr);
FilenameEnd  = sprintf('me_unme_nposm_nposu_meth_End6_%s.txt',name_chr);
FilenameMes  = sprintf('me_unme_nposm_nposu_meth_Mes6_%s.txt',name_chr);

file_countEct=(fullfile(my_folder, FilenameEct));
file_countEnd=(fullfile(my_folder, FilenameEnd));
file_countMes=(fullfile(my_folder, FilenameMes));

[DMR,DMR_ect,DMR_end,DMR_mes,SMRH,SMRL]=compute_dmr_smr_main_2(chrN,name_chr,file_win,file_countEct,file_countEnd,file_countMes,thr_cov,thr_dif,thr_sim,thr_meH,thr_meL);
size(DMR)
size(SMRL)
