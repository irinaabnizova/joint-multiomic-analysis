function [SAR,SARH1,SARL]=chastity_similar_HL_oneL(thr_sim,thr_accH,thr_accL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme);

%---------------do I need all of them
%:SARHe,SARHen,SARHm,SARLe,SARLen,SARLm???

% or just one repr SARH  SARL  ??


   %--------------------------13 fields from 14 May2020
    %  SAR=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' BE_npos_meS' BE_npos_unmeS' BEn_npos_meS' BEn_npos_unmeS' BC_npos_meS' BC_npos_unmeS'];
%           1    2       3       4       5       6        7       8

  %--------------------------13 fields, from 29 June 2019
     % SAR=[sw' sdif' sdif_GC' lev_ES' lev_CS' cov_ES' cov_CS' BE_npos_meS'  BE_npos_unmeS' BC_npos_meS' BC_npos_unmeS' sZ' swi'];

%thr_sim=0.34 for purity

cs=0;
   sw=[];% lowCh/signif win
   swi=[];% lowCh/signif win index
   %sdif=[];
   %sZ=[];
   for i=1:length(indm),
      if pu(i) < thr_sim,
        cs=cs+1;
        sw(cs)=winm(i);% low chastity (equal) wins
        swi(cs)=indm(i);
        %sdif(cs)=dif(i);
        %sdif_GC(cs)=dif_GC(i);
        lev_CS(cs)=levBC(i);
        lev_ES(cs)=levBE(i);
        lev_EnS(cs)=levBEn(i);
        
        cov_CS(cs)=covBC(i);
        cov_ES(cs)=covBE(i);
        cov_EnS(cs)=covBEn(i);
        %sZ(cs)=sigd(i);
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
   
   %--------------------------14 fields
      SAR=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' BE_npos_meS' BE_npos_unmeS' BEn_npos_meS' BEn_npos_unmeS' BC_npos_meS' BC_npos_unmeS' pup'];
%           1    2       3       4       5       6        7       8
%                                                                              9               10            11             12               13         14 
   %--------------------------------SARs high accessible
   
   
      %thr_acc=0.35;

ne=0;nm=0; nen=0;
neen=0;nee=0;nmm=0;% number of low accessed SARs
%SARHm=[]; SARHe=[]; SARHen=[];
%SARLm=[]; SARLe=[];SARLen=[];
SARH=[];SARL=[];



for i=1:length(lev_ES),
    
    %-Ect
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

SARH1=SARH(1:ne-1,:);% changed 3 Sept 2021