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
 
   FilenameEct  = sprintf('me_unme_nposm_nposu_Ect1_%s.txt',name_chr);
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
    
    %================================subfunction prefilter
    function [winn,level_BEFneg,covEFneg,BE_npos_meNeg,BE_npos_unmeNeg,fracZ_lowE]=prefilter_count_chrN_ext(thr_low,FilenameWin,winFolder,my_folder,FilenameEct,chrN)
%filter out low coverage windows for GC methylation
display('%------------------count prefiltering'); 

    win=load(fullfile(winFolder, FilenameWin));
    Y_muEc = load(fullfile(my_folder, FilenameEct));

%  Y_mu is count in a window here, 

%me     umne  #pos_me #pos_unme (5= sum is number of GC)
%0	2	0	1  1
%3	1	1	1  2
%1	1	1	1  2
%0	0	0	0
%0	0	0	0
%0	0	0	0
%2	21	1	4  5
%0	3	0	1
%2	23	0	6
%0	4	0	4
 
    [winn,level_BEFneg,covEFneg,BE_npos_meNeg, BE_npos_unmeNeg,fracZ_lowE]=prefilter_internal_May21(thr_low,win,Y_muEc);

end

%-------------------subfunction internal prefiltering
function [winn,level_BEFneg,covEFneg,BE_npos_meNeg, BE_npos_unmeNeg,fracZ_lowE]=prefilter_internal_May21(thr_low,win,Y_muEc);
%filter out low coverage windows for GC methylation- function itself

   win_len=100;%- hard coded here

%------------------OUTPUTS
%winn,level_BEFneg,covEFneg, BE_npos_meNeg, BE_npos_unmeNeg,
  
%1=================--------choosing slices: for vizualization
    j=1; %just a slice starting from n1=1+N*(j-1), length N=....
    Nint=(max(win)-min(win))/win_len; %1400000 % number of 100 bp windows in the interval
    n1=1+Nint*(j-1);
    n2=Nint*j;%--------------for subregion-10K windows

%2 ========================thr low accessibility FILTERING

display('counts coverage for Layer i');
    [level_BEFneg,covEFneg, BE_npos_meNeg, BE_npos_unmeNeg, level_BEF,covEF,BE_npos_meF,BE_npos_unmeF,winEF, level_BE_low,covE_low,BE_npos_me_low,BE_npos_unme_low,win_LowCovE,winZE,fracZ_lowE]=counts_access_level_nposme_thr_filter_2(thr_low,win,n1,n2,Y_muEc);
    winn=win(n1:n2);

end

%------------------subfunction2: count access
function [level_BMFneg,covMFneg, BM_npos_meNeg, BM_npos_unmeNeg, level_BMF,covMF,BM_npos_meF,BM_npos_unmeF,winMF, level_BM_low,covM_low,BM_npos_me_low,BM_npos_unme_low,win_LowCovM,winZM,fracZ_low]=counts_access_level_nposme_thr_filter_2(thr_low,win,n1,n2,Y_muM);
%  Y is count in a window here, 
    winn=win(n1:n2);
%========================
    BM_me=[Y_muM(n1:n2,1)];%No
    BM_unme=[Y_muM(n1:n2,2)];%Y
    BM_npos_me=[Y_muM(n1:n2,3)];% number/count GC-meth pos
    BM_npos_unme=[Y_muM(n1:n2,4)];% number/count GC-unmeth pos

%------find not covered windows, unmethylated=non-accessible
    BM_mu=Y_muM(n1:n2,:);%[BM_me BM_unme];

%1--------compute level access/methylation---- same here
%2-----------------retrieve uncovered regions

    [level_BM, levCovM,winCovM,BM_meCov, BM_unmeCov,BM_npos_meCov,BM_npos_unmeCov,winZM,propZM]=meth_zero_covered_level_win(BM_me,BM_unme,BM_npos_me,BM_npos_unme,winn);
    [level_BMFneg,covMFneg, BM_npos_meNeg, BM_npos_unmeNeg, level_BMF,covMF,BM_npos_meF,BM_npos_unmeF,winMF, level_BM_low,covM_low,BM_npos_me_low,BM_npos_unme_low,win_LowCovM,prop_lowM]=level_win_low_high_covered_nme_out(winn,level_BM,BM_me,BM_unme,BM_npos_me,BM_npos_unme,thr_low);
    fracZ_low=[propZM,prop_lowM];
end

%-------------------------subfunction
function [level_BE, levCov,winCov,BE_meCov, BE_unmeCov,BE_npos_meCov,BE_npos_unmeCov,winZ,propZ]=meth_zero_covered_level_win(BE_me,BE_unme,BE_npos_me,BE_npos_unme,winn)
%compute level methylation, GC data:  level_BE(i)=BE_me(i)/s(i);% the
%higher - the worse for CpG

%---------------------OUTPUTS
%-------------------covered windows: 
% levCov,winCov,BE_meCov, BE_unmeCov,BE_npos_meCov,BE_npos_unmeCov
%------------------all windows: showing as negatives uncovered windows
%level_BE,   levALL,winALL,BE_meALL, BE_unmeALL,BE_npos_meALL,BE_npos_unmeALL
%1.0-------------separates zero_wins from nonZero: 
%2-----------------retrieve uncovered regions

%--------------------INPUT
%  BE_me,BE_unme  = counts me and unme in each window of the region
%------------winn = windows within the region

%--------------------initialise
    winZ=[];
    winnz=[];
    indz=[];

    winCov=[];
    winnCov=[];
    indCov=[];

    BE_npos_meCov=[];
    BE_npos_unmeCov=[];

%---------------------body

    cz=0;
    c=0;
    for i=1:length(BE_me),
        s(i)=BE_me(i)+BE_unme(i);
    %---------------------nonZ windows
      if s(i)>0,
        c=c+1;
        winnCov(c)=winn(i);
        indCov(c)=i;
        levCov(c)=BE_me(i)/s(i);
        BE_npos_meCov(c)=BE_npos_me(i);
        BE_npos_unmeCov(c)=BE_npos_unme(i);
        BE_meCov(c)=BE_me(i);
        BE_unmeCov(c)=BE_unme(i);
          
        %-------------------all windows
        level_BE(i)=BE_me(i)/s(i);%the higher - the better for chrom access(GC)
        
      else
        cz=cz+1;
        level_BE(i)=-0.1;
        winnz(cz)=winn(i);
        indz(cz)=i;
      end
    end

    cov_uncov_num=[c,cz];
    tot=c+cz;

    winZ=[winnz' indz'];
    winCov=[winnCov' indCov'];

    countZ=cz;
    propZ=cz/tot;
end

%------------------------subfunction level_window...
function [level_BEFneg,covEFneg, BE_npos_meNeg, BE_npos_unmeNeg, level_BEF,covEF,BE_npos_meF,BE_npos_unmeF,winEF, level_BE_low,covE_low,BE_npos_me_low,BE_npos_unme_low,win_LowCov,prop_low]=level_win_low_high_covered_nme_out(winn,level_BE,BE_me,BE_unme,BE_npos_me,BE_npos_unme,thr_low)

%======filter out low and high number 
%---------------------------------------addded 
%BE_npos_me,BE_npos_unme  = for each window, number of GCmeth and GC-unmeth
%positions (max 50)

%-----------------------------OUTPUTs
% all about filtered inIN: level, cov,num_meth, num_unmeth,windows
%level_BEF,covEF,BE_npos_meF,BE_npos_unmeF,winEF,
% all about filtered OUT: level, cov,num_meth, num_unmeth,windows
%level_BE_low,covE_low,BE_npos_me_low,BE_npos_unme_low,win_LowCov,prop_low
%-----------------------output with deliberate negatives
  %level_BEFneg(i)=-0.5;%level_BE(i);
   %     covEFneg(i)=-300;%covE(i);
   
   %   level_BEF(cfi)=level_BE(i);-------Filtered: without empty and low-covered  win
   %     covEF(cfi)=covE(i);---------coverage
   
%---------------------------------initialise
    ind_low=[];
    win_low=[];
    indF=[];winF=[];
    level_BEF=[];covEF=[];
    level_BE_low=[];covE_low=[];
%------------------------------------body
    cfi=0;% count filtered in
    cfo=0;
    for i=1:length(level_BE),
        covE(i)=BE_me(i)+BE_unme(i);
        if covE(i)>=thr_low, %leave it
            cfi=cfi+1;
            level_BEF(cfi)=level_BE(i);
            covEF(cfi)=covE(i);
            indF(cfi)=i;
            winF(cfi)=winn(i);
            BE_npos_meF(cfi)=BE_npos_me(i); 
            BE_npos_unmeF(cfi)=BE_npos_unme(i); 
            %---------for plotting , with all (negative including) windows
            level_BEFneg(i)=level_BE(i);
            covEFneg(i)=covE(i);
            BE_npos_meNeg(i)=BE_npos_me(i); 
            BE_npos_unmeNeg(i)=BE_npos_unme(i); 
        else
            cfo=cfo+1;
            ind_low(cfo)=i;
            win_low(cfo)=winn(i);
        
            level_BEFneg(i)=-0.5;%level_BE(i);
            covEFneg(i)=-300;%covE(i);
            BE_npos_meNeg(i)=-150; 
            BE_npos_unmeNeg(i)=-150; 
        
            BE_npos_me_low(cfo)=BE_npos_me(i); 
            BE_npos_unme_low(cfo)=BE_npos_unme(i); 
            level_BE_low(cfo)=level_BE(i);
            covE_low(cfo)=covE(i);
        end
     end

    winEF=[indF' winF'];
    win_LowCov=[ind_low' win_low'];

    count_low=i-cfi;
    prop_low=(i-cfi)/i;
    highCov_lowCov=[cfi,cfo];

end
   
%=================================subfunction save DAR
    function [chrN]=save_DAR_EctEndMes(DAR,folder,textFilename,chrN)
%===============================write DAR or SAR in the 8 field's file
%------------------------Input
 %SAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %    1    2  3  4     5        6     7       8
    chr=DAR(:,1);
    starts=DAR(:,2)-25;
    ends=DAR(:,3)+25;
    ch=DAR(:,4);
    lev_EcS=DAR(:,5);
    lev_EnS=DAR(:,6);
    lev_MS=DAR(:,7);
    ind=DAR(:,8);
       
  %--------------------DARs 8 fileds
  %textFilename=sprintf('SARL_EctEndMes_c25_aL02_GCPu.txt');
  fp = fopen(fullfile(folder, textFilename), 'w');
  for i=1:length(starts),
       fprintf(fp,'%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n',chr(i),starts(i),ends(i),ch(i),lev_EcS(i),lev_EnS(i),lev_MS(i),ind(i));
  end
  fclose(fp);
 end

