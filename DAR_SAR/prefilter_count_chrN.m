function [winn,level_BEFneg,covEFneg,BE_npos_meNeg,BE_npos_unmeNeg,fracZ_lowE]=prefilter_count_chrN(thr_low,FilenameWin,winFolder,my_folder,FilenameEct)
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

%=====================------subfunction: count access
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

%===================--subfunction: filter Zero covered windows
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

%================--------subfunction: preprocess level_windows (mb for
%visual as well)
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
