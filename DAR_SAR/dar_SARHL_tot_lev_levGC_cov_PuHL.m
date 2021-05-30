function [tot, DAR, DAR_raw, SAR_raw,SARH,pu,ch,amp,SARL]=dar_SARHL_tot_lev_levGC_cov_PuHL(thr_sim,thr_dif,thr_accH,thr_accL,winn,level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg, chrN);

% computes DARs SARs due to Chastity/Purity in one go
% SARHL- separate H L accessibility thr: thr_accH,thr_accL

%------------------OUTPUTs important are DAR, SARH:format
 %SARH=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %    1    2  3  4     5        6     7       8

%----------------- INPUTS: precomputed FILTERED counts of GC me unme in a 100bp
%BC_npos_unmeNeg,BE_npos_meNeg, BC_npos_unmeNeg  precomputed from the :
%[level_BMFneg,covMFneg, BM_npos_meNeg, BM_npos_unmeNeg, level_BMF,covMF,BM_npos_meF,BM_npos_unmeF,winMF, level_BM_low,covM_low,BM_npos_me_low,BM_npos_unme_low,win_LowCovM,winZM,fracZ_lowM]=counts_access_level_nposme_thr_filter_1(thr_low,win,n1,n2,Y_muM);

%----------------------Params
    thr_sim_dif=[thr_sim,thr_dif]% for Pu
    thrHL_access=[thr_accH,thr_accL]% for accessibility level: high is good/active

%-------------initialise
    tot=[]; % all mutually decently covered
    pu=[];ch=[];
    amp=[];
    DAR=[];DAR_raw=[];
    SARH=[];SARL=[];

  
%1------------find mutually non-empty decently covered windows
    %winn=win(n1:n2);
    
    cm=0;% count mutual
    for i=1:length(level_BEFneg),
        
        if level_BEFneg(i)>=0 && level_BCFneg(i)>=0 && level_BEnFneg(i)>=0,% if ==0, it is low covered of zero
            cm=cm+1;
            indm(cm)=i;
            winm(cm)=winn(i);
            levBE(cm)=level_BEFneg(i);
            levBC(cm)=level_BCFneg(i);
            levBEn(cm)=level_BEnFneg(i);
            
            covBE(cm)=covEFneg(i);
            covBC(cm)=covCFneg(i);
            covBEn(cm)=covEnFneg(i);
  
            BE_npos_me(cm)=BE_npos_meNeg(i);            
            BC_npos_me(cm)=BC_npos_meNeg(i);
            BEn_npos_me(cm)=BEn_npos_meNeg(i);
             
            BE_npos_unme(cm)=BE_npos_unmeNeg(i);
            BC_npos_unme(cm)=BC_npos_unmeNeg(i);
            BEn_npos_unme(cm)=BEn_npos_unmeNeg(i);
            
            levGC_BE(cm)=BE_npos_me(cm)/(BE_npos_me(cm)+BE_npos_unme(cm));
            levGC_BC(cm)=BC_npos_me(cm)/(BC_npos_me(cm)+BC_npos_unme(cm));
            levGC_BEn(cm)=BEn_npos_me(cm)/(BEn_npos_me(cm)+BEn_npos_unme(cm));
      
        end
    end
    
      
 %2===================Stats Levels and  GC across windows
    
    if cm >0,% continue if any wins are selected
   
 %2.1======================define total for mutually OK windows
 
     %-----------------------------------------14(with GC) fieleds
     tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me'  BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%    tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%        1     2      3      4        5      6       7       8          9
%                                                                                      10          11              12            13     14 
 %3=====================compute Chastity/Purity of tot(:,2:4) (three
     %access Levels
    
     [pu,ch,amp]=get_pu_ch_access(tot);
     
 %4  DARs ========================find dominant windows via Pu Ch, thr_dif
     [DAR,DAR_raw]=chastity_differential_format(thr_dif,thr_accH,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme,chrN);
 
 %5=======================================similars: HL is di-chotomus
   
     [SAR,SARH,SARL,SAR_raw]=chastity_similar_format(thr_sim,thr_accH,thr_accL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme,chrN);
 
   end % if cm>0 any mutually covered windows exist
end % function

%-------------------------subfunctions 
function [pu,ch,amp]=get_pu_ch_access(tot_access)

%compute Purity and Chastity for each non empty access window
%----------------desinged for three layers currently, e.g. E7.5 ect end mes

%============================INPUTS
%tot_access=   
%tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me'  BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%  tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%        1     2      3      4        5      6       7       8          9
%                                                                                      10          11              12            13     14 
    pu=[];
    ch=[];
    amp=[];

%1----------------------the first round, to learn real amplitude
    layers=tot_access(:,2:4);% levels of accessibility
    nn=0;
    for i=1:length(tot_access(:,1)),
        vec=layers(i,:);
        mi=min(vec);
        if mi > 0,
        pu(i)=max(vec)/sum(vec);
        amp(i)=max(vec);
        ch(i)=max(vec)/(max(vec)+second_max(vec));
        else
        nn=nn+1;
        con=abs(mi);
        amp(i)=max(vec);
        vec1=vec+con;
        if sum(vec1)>0,
        pu(i)=max(vec1)/sum(vec1);
        ch(i)=max(vec1)/(max(vec1)+second_max(vec1));
        else
        pu(i)=0;
        ch(i)=0;
        end
    end
end
    

figure;
subplot(2,1,1);
histogram(pu);
title('pu across three layers');
xlabel('Pu of GE across Ect End Mes, each access win');
xlim([0 1]);
subplot(2,1,2);
histogram(ch);
title('chastity across three layers');
xlim([0 1]);
 
end

%----------------------subfunction second max
function [second_max]=second_max(vecc)
   ma=max(vecc);
   vecc1=vecc;
   for i=1:length(vecc)
      if  vecc1(i)==ma,
        vecc1(i)=min(vecc);
      end
   end

   second_max=max(vecc1);
end

%---------------------subfunction get dars 
function [DAR,DAR_raw]=chastity_differential_format(thr_dif,thr_accH,ch,pu,amp,winm,indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme,chrN);

%for all windows computes dominance index for accessibility level
% retrieves windows where one lineage dominates two oters (high
% Pu/Chastity)== DARs
%-----------------done for one chromosome only

%------------------------OUTPUT
 %DAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %    1    2  3  4     5        6     7       8

%---------------DARs for any of three layers
    win_len=100;%hard coded 

   cs=0;
   sw=[];% lowCh/signif win
   swi=[];% lowCh/signif win index
   for i=1:length(indm),
      if pu(i) > thr_dif & amp(i)>=thr_accH,% max acc sh be large enough
        cs=cs+1;
        sw(cs)=winm(i);% high chastity (dominant) wins
        swi(cs)=indm(i);
        lev_CS(cs)=levBC(i);
        lev_ES(cs)=levBE(i);
        lev_EnS(cs)=levBEn(i);
        
        cov_CS(cs)=covBC(i);
        cov_ES(cs)=covBE(i);
        cov_EnS(cs)=covBEn(i);
        
            BE_npos_meS(cs)=BE_npos_me(i);            
            BC_npos_meS(cs)=BC_npos_me(i);
            BEn_npos_meS(cs)=BEn_npos_me(i);     
            
            BE_npos_unmeS(cs)=BE_npos_unme(i);
            BC_npos_unmeS(cs)=BC_npos_unme(i);
            BEn_npos_unmeS(cs)=BEn_npos_unme(i);
            chp(cs)=ch(i);
            pup(cs)=pu(i);
 
      end
   end
   
%---------------sum up GC counts
   numGC_ec=BE_npos_meS' + BE_npos_unmeS';
   numGC_en=BEn_npos_meS' + BEn_npos_unmeS';
   numGC_m=BC_npos_meS' + BC_npos_unmeS';
   ind=1:length(chp);
  
%      DAR_raw=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' BE_npos_meS' BE_npos_unmeS' BEn_npos_meS' BEn_npos_unmeS' BC_npos_meS' BC_npos_unmeS' chp' pup'];
%           1    2       3       4       5       6        7       8
%                                                                              9               10            11             12               13           14     15
  
    DAR_raw=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' numGC_ec numGC_en numGC_m chp' pup' ind'];
   %          1    2       3       4       5       6        7       8
   %                                                                           9        10     11  12   13       

%===================make more standard format
   
   chr=chrN*ones(length(chp),1);% column
   st=sw';
   en=sw'+ win_len;
   DAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %    1    2  3  4     5        6     7       8
end

%------------------------subfunction get SARs
function [SAR,SARH,SARL,SAR_raw]=chastity_similar_format(thr_sim,thr_accH,thr_accL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme,chrN);
%for all windows computes dominance index for accessibility level
% retrieves windows where one lineage similar to two others (low
% Pu/Chastity)== SARs
%-----------------done for one chromosome only

%------------------------OUTPUT
 %SAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %    1    2  3  4     5        6     7       8
   win_len=100;%hard coded
%thr_sim=0.34 for low purity

   cs=0;
   sw=[];% lowCh/signif win
   swi=[];% lowCh/signif win index
   for i=1:length(indm),
      if pu(i) < thr_sim,
        cs=cs+1;
        sw(cs)=winm(i);% low chastity (equal) wins
        swi(cs)=indm(i);
        lev_CS(cs)=levBC(i);
        lev_ES(cs)=levBE(i);
        lev_EnS(cs)=levBEn(i);
        
        cov_CS(cs)=covBC(i);
        cov_ES(cs)=covBE(i);
        cov_EnS(cs)=covBEn(i);
        %-------------------------GC counts
            BE_npos_meS(cs)=BE_npos_me(i);            
            BC_npos_meS(cs)=BC_npos_me(i);
            BEn_npos_meS(cs)=BEn_npos_me(i);     
            
            BE_npos_unmeS(cs)=BE_npos_unme(i);
            BC_npos_unmeS(cs)=BC_npos_unme(i);
            BEn_npos_unmeS(cs)=BEn_npos_unme(i);
            chp(cs)=ch(i);
            pup(cs)=pu(i);
 
      end
   end
   
  %--------------------------sum up GCs
   numGC_ec=BE_npos_meS' + BE_npos_unmeS';
   numGC_en=BEn_npos_meS' + BEn_npos_unmeS';
   numGC_m=BC_npos_meS' + BC_npos_unmeS';
   ind=1:length(chp);
  
%      DAR_raw=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' BE_npos_meS' BE_npos_unmeS' BEn_npos_meS' BEn_npos_unmeS' BC_npos_meS' BC_npos_unmeS' chp' pup'];
%           1    2       3       4       5       6        7       8
%                                                                              9               10            11             12               13           14     15
  
      SAR_raw=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' numGC_ec numGC_en numGC_m chp' pup' ind'];
   %           1    2       3       4       5       6        7       8
   %                                                                           9        10     11  12   13       
%===================make more standard format
      chr=chrN*ones(length(chp),1);% column
      st=sw';
      en=sw'+ win_len;
      SAR=[chr st en chp' lev_ES' lev_EnS' lev_CS' ind'];
    %       1    2  3  4     5        6     7       8

%--------------------------------SARs high accessible
         %thr_acc=0.35;
ne=0;
nee=0;% number of low accessed SARs
SARH=[];SARL=[];

for i=1:length(lev_ES),
    
    if lev_ES(i)>= thr_accH && lev_EnS(i)>= thr_accH && lev_CS(i)>= thr_accH,
        ne=ne+1;
        SARH(ne,:)=SAR(i,:);
    end
    if lev_ES(i)<= thr_accL && lev_EnS(i)<= thr_accL && lev_CS(i)<= thr_accL,
        nee=nee+1;
        SARL(nee,:)=SAR(i,:);
    end
end
nHL=[ne,nee]

end
   