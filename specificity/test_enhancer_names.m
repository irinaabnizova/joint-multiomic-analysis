function [FLC,list_SEG,list_Ect,list_End,list_Mes,outOriginal,outShuffled]=test_enhancer_names()
   
  
   my_folder='data_specificity\';
   folder='outputs\';
   
   %FLC= dictionary of all motifs
   
%cd data_specificity\
%ls

%                              dmr_Mes_1.m                    
%..                             feature_data_num.txt           
%SARH_1.m                       normal_feature_name_long_2.m   
%capital_feature_name_long_2.m  p_SEG_1.m                      
%dar_ect.m                      p_ect_1.m                      
%dar_end.m                      p_end_1.m                      
%dar_mes.m                      p_mes_1.m                      
%dmr_Ect_1.m                    smrL_SEG_1.m                   
%dmr_End_1.m                    

    
   FilenameSEG=sprintf('SARH_1.m');
   FilenameEct  = sprintf('dar_ect.m');
   FilenameEnd  = sprintf('dar_end.m');
   FilenameMes  = sprintf('dar_mes.m');
   
   nameMatrixOriginal=sprintf('Matrix_SARH_EctEndMes_output.txt');
   nameShuffled=sprintf('Matrix_shuffled_output.txt');
     
   
        nameDictionary=sprintf('capital_feature_name_long_2.m');

       list_SEG=fullfile(my_folder, FilenameSEG);
       list_Ect=fullfile(my_folder, FilenameEct);
       list_End=fullfile(my_folder, FilenameEnd);
       list_Mes=fullfile(my_folder, FilenameMes);
       
       FLC=fullfile(my_folder,nameDictionary);
       
       outOriginal=fullfile(folder,nameMatrixOriginal);
       outShuffled=fullfile(folder, nameShuffled);
       
       
       