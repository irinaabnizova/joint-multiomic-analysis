function [file_win,file_countEct,file_countEnd,file_countMes,outDAR,outSAR]=test_file_names(name_chr)
   
   winFolder=('test_data\windows\');
   my_folder='test_data\counts\';
   folder='test_data\output_DARSAR\';
    
   FilenameWin=sprintf('wine100_%s.txt',name_chr);
   FilenameEct  = sprintf('me_unme_nposm_nposu_Ect1_%s.txt',name_chr);
   FilenameEnd  = sprintf('me_unme_nposm_nposu_End5_%s.txt',name_chr);
   FilenameMes  = sprintf('me_unme_nposm_nposu_Mes_%s.txt',name_chr);
   
   nameSAR=sprintf('SARH_EctEndMes_output_%s.txt',name_chr);
   nameDAR=sprintf('DAR_EctEndMes_output_%s.txt',name_chr);
 

       file_win=fullfile(winFolder, FilenameWin);
       file_countEct=fullfile(my_folder, FilenameEct);
       file_countEnd=fullfile(my_folder, FilenameEnd);
       file_countMes=fullfile(my_folder, FilenameMes);
       outSAR=fullfile(folder,nameSAR);
       outDAR=fullfile(folder,nameDAR);
       