function [DMR,DMR_ect,DMR_end,DMR_mes]=just_DMR_thrSept(thr_dif,thr_meH,thr_meL,ch,pu,amp,winm, indm,levBE,levBEn,levBC,covBE,covBEn,covBC,BE_npos_me,BEn_npos_me,BC_npos_me,BE_npos_unme,BEn_npos_unme,BC_npos_unme);

%-Sept 2021- to increase DMR of Ect mes (because End though less abindand
%non empty windows, but lower meth within them!!!)

%thr_accH=thr_meH
%thr_accL=thr_meL+0.4
%thr_dif=0.35 % for ch/pu

%thr_sim=0.34 for purity

%-------------------------outputs
 %  DMR=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' pu'];
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

