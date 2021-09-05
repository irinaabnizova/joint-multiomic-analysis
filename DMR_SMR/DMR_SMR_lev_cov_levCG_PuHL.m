function [tot1,SMR,SMRH,SMRL,DMR,DMR_ect,DMR_end,DMR_mes,pu,ch,amp]=DMR_SMR_lev_cov_levCG_PuHL(thr_dif,thr_sim,thr_meH,thr_meL,winn,level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg);

% Sept 2021 changed just_DMR_thr into just_DMR_thrSept (%-Sept 2021- to increase DMR of Ect mes (because End though less abindand
%non empty windows, but lower meth within them!!!)


%-----------------for sorted wins in one chromosome
%------------------------------------due to Pu (highest dominance 1-me)
%                                                                                     10             11           12           13          14                                   
%SARHe,SARLe  --- thrHL expected will be more or
%less the same for any Layer!!!...

% DAR SAR due to Chastity/Purity in one go
% SARHL- separate H L access thr: thr_meH,thr_meL

%-----------------------INPUTS

%   thr_meH,thr_meL,
%thr_dif,thr_sim,thr_meH,thr_meL,

%winn,
%----------------------levels of DNA meth and numbers of me-unme positions
%level_BEFneg,level_BEnFneg,level_BCFneg,covEFneg,covEnFneg,covCFneg,BE_npos_meNeg,BEn_npos_meNeg,BC_npos_meNeg, BE_npos_unmeNeg,BEn_npos_unmeNeg,BC_npos_unmeNeg

%----------------------OUTPUTs

%----------------------------output tot1, DMR, SMR 
 %-----------------------------------------14 fields
%  tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%        1     2      3      4        5      6       7       8          9



%=====================specifically for DMRs- even less contrast
thr_meH1=thr_meH;%0.6
thr_meL1=thr_meL+0.15;

thr_SD_pu=[thr_sim,thr_dif]
%thr_sim=0.51
thrHL_meS=[thr_meH,thr_meL]
thrHL_meD=[thr_meH1,thr_meL1]% less contrast than Sim

 % DAR_car=[starts ends lev_ES' cov_ES' numGCe lev_CS' cov_CS' numGCm sdif'
 % ind'];no chr here

  %-INPUT  =  level_BEFneg, level_BCFneg - 
  %  Prepresecced: if a window poorely covered, it is declared negative value
  
  tot=[]; % all mutually decently covered across three layers
  pu=[];
  ch=[];
  amp=[];
  DMR=[]; sw=[];% all signif 
  SMR=[];
  %DMR_car=[];
  SMRH=[];SMRL=[];

  
 %1===========-------find mutually non-empty decently covered windows
    %winn=win(n1:n2);
    
    cm=0;% count mutual windows
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
             
            %BE_npos_meNeg, BC_npos_meNeg, BE_npos_unmeNeg,BC_npos_unmeNeg
            BE_npos_me(cm)=BE_npos_meNeg(i);            
            BC_npos_me(cm)=BC_npos_meNeg(i);
            BEn_npos_me(cm)=BEn_npos_meNeg(i);
             
            BE_npos_unme(cm)=BE_npos_unmeNeg(i);
            BC_npos_unme(cm)=BC_npos_unmeNeg(i);
            BEn_npos_unme(cm)=BEn_npos_unmeNeg(i);
            
            levGC_BE(cm)=BE_npos_me(cm)/(BE_npos_me(cm)+BE_npos_unme(cm));
            levGC_BC(cm)=BC_npos_me(cm)/(BC_npos_me(cm)+BC_npos_unme(cm));
            levGC_BEn(cm)=BEn_npos_me(cm)/(BEn_npos_me(cm)+BEn_npos_unme(cm));
            
     end % if mutially covered
 end
    
      
    %===================
 if cm >0,% continue if any wins are selected
           
  
   %2===============================define total for mutually OK windows
      %-----------------------------------------14(with GC) fieleds
     tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me'  BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
 %  tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%        1     2      3      4        5      6       7       8          9
%                                                                                      10          11              12            13     14 
    %3=====================compute Chastity/Purity of tot(:,2:4) (three
    %meth Levels
    %----------------compute pu of me1=1-me good opened if high (to
    %apply Pu- dominance index)
    [pu,ch,amp]=get_pu_ch_meth_3(tot);
     
    tot1=[tot pu']; % one more field n15
                                                                                       
      %4=======================================similars: HL is di-chotomus
      [SMR,SMRH,SMRL]=chastity_similar_HL_oneL(thr_sim,thr_meH,thr_meL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme);

      %5=============differentially methylated: less strict difference than
      %similarity
   
      %--------------------------14 fields outputs
     % DMR=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' BE_npos_meS' BE_npos_unmeS' BEn_npos_meS' BEn_npos_unmeS' BC_npos_meS' BC_npos_unmeS' pup'];
%           1    2       3       4       5       6        7       8
%    
      
%5.1---there is the endoderm-specific adjustment of HL meth thr further
      [DMR,DMR_ect,DMR_end,DMR_mes]=just_DMR_thrSept(thr_dif,thr_meH1,thr_meL1,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme);

%
  end % if cnm
end
%2 subfuction

  %=====================================SUBFUNCTIONS forr DMR_SMR_lev_cov_levCG_PuHL
function [DMR,DMR_ect,DMR_end,DMR_mes]=just_DMR_thrSept(thr_dif,thr_meH,thr_meL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme);

%-Sept 2021- to increase DMR of Ect mes (because End though less abindand
%non empty windows, but lower meth within them!!!)

%thr_accH=thr_meH
%thr_accL=thr_meL+0.4
%thr_dif=0.35 % for ch/pu

%thr_sim=0.34 for purity

%-------------------------outputs
 %  DMR=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%        1     2      3      4        5      6       7       8          9
%                                                                                      10          11              12            13     14 


%0===================initiate
DMR=[];DMR_ect=[];DMR_end=[];DMR_mes=[];
   cs=0;% count signif=dominant
   sw=[];% lowCh/signif win
   swi=[];% lowCh/signif win index
   lev_CS=[];
   lev_ES=[];
   lev_EnS=[];
    
 %1=======================select windows with high Pu (dominant)  
   for i=1:length(indm),
      if pu(i) > thr_dif,
           cs=cs+1;
           sw(cs)=winm(i);% high pu wins
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
   
   %--------------------------14 fields
      DMR=[sw' lev_ES' lev_EnS' lev_CS' cov_ES' cov_EnS' cov_CS' BE_npos_meS' BE_npos_unmeS' BEn_npos_meS' BEn_npos_unmeS' BC_npos_meS' BC_npos_unmeS' pup'];
%           1    2       3       4       5       6        7       8
%                                                                              9               10            11             12               13        14   
   
  %3=========================sort them by layers (automatically filyering!)
    thr_meH1=thr_meH-0.05;%=0.6500 % more for ect mes
    thr_meL1=thr_meL-0.1;%+0.4500 % less for end

    DMR_ect=[];
    DMR_end=[];
    DMR_mes=[];
    ne=0;nee=0;nc=0;
   for i=1:length(lev_ES),
    
        %-------------------DMR Ect: low meth in Ect, high meth in End Mes
        if lev_ES(i)<= thr_meL && lev_EnS(i)> thr_meH1 && lev_CS(i)> thr_meH1,
            ne=ne+1;
            DMR_ect(ne,:)=DMR(i,:);
        end
    
        %------------------------DMR End (more outliers because less
    %coverage?....
       if lev_EnS(i)<= thr_meL1 && lev_ES(i)> thr_meH && lev_CS(i)> thr_meH,
           nee=nee+1;
           DMR_end(nee,:)=DMR(i,:);
       end
    
      %------------------------DMR Mes
      if lev_CS(i)<= thr_meL && lev_EnS(i)> thr_meH1 && lev_ES(i)> thr_meH1,
          nc=nc+1;
          DMR_mes(nc,:)=DMR(i,:);
      end
    
   end
   numDMR=[ne,nee,nc]

end  
 %============================subfunction 2 
function [SAR,SARH1,SARL]=chastity_similar_HL_oneL(thr_sim,thr_accH,thr_accL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme);

% or just one repr SARH  SARL  
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
end