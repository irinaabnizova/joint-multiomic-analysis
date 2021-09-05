function [chrN]=save_DMR_SMR_EctEndMes_H07_L05(DMR_ect,DMR_end,DMR_mes,SMRH,SMRL,name_chr,chrN)
%===============================save DARi SARH SARL without GC agreement

%---------------------------INPUT tot1, DMR, SMR June 2020
 %-----------------------------------------14(with CG and pu) fieleds
%  tot1=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' pu'];
%        1     2      3      4        5      6       7       8          9
%        d     .2f    .2f    .2f      d      d       d       .2f          d               10             11            12            13         14   
%                                                                                        d             d              d             d        .2f                                                                      
%filtering
chrN

%====================what format? Carines but no flanks?
% % 10 fields ( 10 with chr)
  %DAR_car=[chr starts ends lev_ES  cov_ES numGCe lev_CS cov_CS numGCm dif ind];


  %folder='C:\Users\ia1\Documents\MATLAB\DMR_SMR_OPSB\data_dmr_opsb\dmr_smr_unique\';
  folder=('test_data\output_DMRSMR\');
 
  %[SAR_car]=format_SAR_carine(flank,SAR, chrN);
  %[DAR_car]=format_DAR_carine(flank,DAR, chrN);
  
   
  %--------INPUT sets in carines format
    % 11 fields ( 11 with chr)
  %DAR_car=[chr starts ends lev_ES  cov_ES numGCe lev_CS cov_CS numGCm dif ind];
  %          1   2      3     4      5      6       7      8     9      10
  %          11
     %    13.00    3713498.00    3713623.00          0.57         49.00          5.00          0.20
     %    13.00    3878098.00    3878223.00          0.40         30.00          6.00          0.07
       
     %      30.00          5.00          0.37         28.00
    %     28.00          6.00          0.33         42.00
    
    
    
    
    %========================NOW
     
    %[   chr start    end  levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' pu'];
%        1     2      3      4        5      6       7       8          9         10            11              12            13         14            15

%        d     d      d     .2f     .2f    .2f       d       d          d          d             d              d             d          d             .2f

  %--------------------DARect CG 16 fields with chrN st en
  textFilename=sprintf('DMR_Ect_c4_H07_L05_pu05_%s.txt',name_chr);
  fp = fopen(fullfile(folder, textFilename), 'w');

for i=1:length(DMR_ect(:,1)),
           fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',chrN,DMR_ect(i,1),DMR_ect(i,1)+200,DMR_ect(i,2),DMR_ect(i,3),DMR_ect(i,4),DMR_ect(i,5),DMR_ect(i,6),DMR_ect(i,7),DMR_ect(i,8),DMR_ect(i,9),DMR_ect(i,10),DMR_ect(i,11),DMR_ect(i,12),DMR_ect(i,13),DMR_ect(i,14));
end
 fclose(fp);
 
 %----------------End DMR
 
  textFilename=sprintf('DMR_End_c3_H07_L05_pu05_%s.txt',name_chr);
  fp = fopen(fullfile(folder, textFilename), 'w');

for i=1:length(DMR_end(:,1)),
      fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',chrN,DMR_end(i,1),DMR_end(i,1)+200,DMR_end(i,2),DMR_end(i,3),DMR_end(i,4),DMR_end(i,5),DMR_end(i,6),DMR_end(i,7),DMR_end(i,8),DMR_end(i,9),DMR_end(i,10),DMR_end(i,11),DMR_end(i,12),DMR_end(i,13),DMR_end(i,14));
end
 fclose(fp);
 
 
 %----------------Mes DMR
 
  textFilename=sprintf('DMR_Mes_c4_H07_L05_pu05_%s.txt',name_chr);
  fp = fopen(fullfile(folder, textFilename), 'w');

for i=1:length(DMR_mes(:,1)),
       fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',chrN,DMR_mes(i,1),DMR_mes(i,1)+200,DMR_mes(i,2),DMR_mes(i,3),DMR_mes(i,4),DMR_mes(i,5),DMR_mes(i,6),DMR_mes(i,7),DMR_mes(i,8),DMR_mes(i,9),DMR_mes(i,10),DMR_mes(i,11),DMR_mes(i,12),DMR_mes(i,13),DMR_mes(i,14));
end
 fclose(fp);
 
 
 %----------------SMRH- all methylated! >70%
 
  textFilename=sprintf('SMRH_c4_H07_pu05_%s.txt',name_chr);
  fp = fopen(fullfile(folder, textFilename), 'w');

for i=1:length(SMRH(:,1)),
       fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',chrN,SMRH(i,1),SMRH(i,1)+200,SMRH(i,2),SMRH(i,3),SMRH(i,4),SMRH(i,5),SMRH(i,6),SMRH(i,7),SMRH(i,8),SMRH(i,9),SMRH(i,10),SMRH(i,11),SMRH(i,12),SMRH(i,13),SMRH(i,14));
end
 fclose(fp);
 
 %----------------SMRL un-meth <50%
 
  textFilename=sprintf('SMRL_c4_L05_pu05_%s.txt',name_chr);
  fp = fopen(fullfile(folder, textFilename), 'w');

for i=1:length(SMRL(:,1)),
       fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',chrN,SMRL(i,1),SMRL(i,1)+200,SMRL(i,2),SMRL(i,3),SMRL(i,4),SMRL(i,5),SMRL(i,6),SMRL(i,7),SMRL(i,8),SMRL(i,9),SMRL(i,10),SMRL(i,11),SMRL(i,12),SMRL(i,13),SMRL(i,14));
end
 fclose(fp);
 
 