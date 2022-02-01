function [ssm2,CAC_opt,ng_LMH,ne_LMH,prop_LMH,cGE1,cGE2,cGE3,tpp,tppg,CACC,dCC,den_lin]=get_opt_CAC_den_perf_chrN(thr_vic_array,thrh,t_GE,te,gene_names,ectN,endN,mesN,indEctN,indEndN,indMesN,vec_EcN,vec_EnN,vec_MN,dar1,dar2,dar3,chrN);

%  was part of m1_m2_50K_chrN_11
%---------------------so far outputs only opt CAC: to add second pass:
%                    second pass: get linked enh and genes according to this CAC opt 
%------------------------------make a loop on k: maximise density and cAC

 %  m2_enh=ma_ect+ma_end+2*ma_mes;%????
 CAC_opt=0;
 
 %=======================params
   thr10=0;
   factor=2000;% to compute number of enhancers, each 2000 bp long


%------------outputs
  %ssm1=sens_spec_m1_chrN;
  %ssm2=sens_spec_m2_chrN;
   ssm2=[0,0];
   
   
   
   %============================body

      n=length(thr_vic_array);
   for k=1:n
          %thr10=0;
          thr_vic=thr_vic_array(k);
          
          [ng_LMH,ne_LMH,prop_LMH,TabH1,TabH2,TabH3,cGE1,cGE2,cGE3]=num_enh_lineage_11(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,endN,mesN,indEctN,indEndN,indMesN,vec_EcN,vec_EnN,vec_MN,dar1,dar2,dar3,chrN);
          %prop_LMH
 %--CAC-based decision1--------we accept per lineage in case last row is more than first row
%prop_LMH =

 %         1.44          1.09          0.41
 %         1.33          1.20          0.33
 %         1.67          1.12          1.40
 
         [CAC,tp,fp,TPe,FPe,FNe,tpg,fpg,TPg,FPg,FNg]=CAC_simplistic_decisionLineage(prop_LMH,ng_LMH,ne_LMH);
 
        % display([' M2: The first CAC decision  TP and FP FN , at this vicinity=',num2str(thr_vic)]);
    
        %===========================optimise dC=den*CAC wihtin lineage
        %j=1,2,3 ect end mes
        for j=1:length(CAC),
            den_lineage(j)=tp(j)/(thr_vic/factor);
            dC(j)=den_lineage(j)*CAC(j);
        end
        %dC
        
        %=across three lineages, for vic(k)
        density_M2=TPe/(thr_vic/factor);
        TNe=thr_vic/1000 -(TPe+FPe); % can be negative...
        perform_m2=[TPe,FPe,TNe,FNe];
     
       %
        DEN(k)=density_M2;
        perf(k,:)=perform_m2;
        
        %====================per lineage for vc(k)
        CACC(k,:)=CAC;% per lineage
        tpp(k,:)=tp;
        fpp(k,:)=fp;
        tppg(k,:)=tpg;
        fppg(k,:)=fpg;
     
        den_lin(k,:)=den_lineage;
        dCC(k,:)=dC;%optimal density*CAC
        
        
        prop_LMHH(:,:,k)=prop_LMH;
        ne_LMHH(:,:,k)=ne_LMH;
        ng_LMHH(:,:,k)=ng_LMH;
        %m1(k)=sum(ne_LMH(3,:));
         
   end % k per vic(k)
        
    %========we optimise min vi, first max TP (found Enh when CAC is
    %posive)
     
%tpp =
% enh ect  end   mes per vic   
  %   0    36     8
  %   0    45    10
  %   0    52    11
  %   0    64    11
  %   7    68    11
  
  %=====optimal within each lineage:
  tp_ect=tpp(:,1);
  CACC_ect=CACC(:,1);
  dCC_ect=dCC(:,1);
  den_lin_ect=den_lin(:,1);
  [ zone_ect,ma_tp_ect,ma_dC_ect,CAC_opt_ect,ind_vic_ect] = opt_vic_CAC( tp_ect,thr_vic_array,CACC_ect,dCC_ect,den_lin_ect);
%  [ vic_ect,ma_ect ] = opt_vic( tp_ect,thr_vic_array );
  
  tp_end=tpp(:,2);
  CACC_end=CACC(:,2);
  dCC_end=dCC(:,2);
  den_lin_end=den_lin(:,2);
  [ zone_end,ma_tp_end,ma_dC_end,CAC_opt_end,ind_vic_end] = opt_vic_CAC( tp_end,thr_vic_array,CACC_end,dCC_end,den_lin_end);
%  [ vic_end,ma_end ] = opt_vic( tp_end,thr_vic_array );
  
   tp_mes=tpp(:,3);
  %[ vic_mes,ma_mes ] = opt_vic( tp_mes,thr_vic_array );
   CACC_mes=CACC(:,3);
   dCC_mes=dCC(:,3);
   den_lin_mes=den_lin(:,3);
   [ zone_mes,ma_tp_mes,ma_dC_mes,CAC_opt_mes,ind_vic_mes] = opt_vic_CAC( tp_mes,thr_vic_array,CACC_mes,dCC_mes,den_lin_mes);

   %====================compute combined performance
    
  fa=1;
  m2_enh=ma_tp_ect+ma_tp_end+fa*ma_tp_mes;%????
  vic_m2=round(mean([zone_ect,zone_end,zone_mes]));
  
  %----------------------corresponds to thr_vic=50K
  %vic_m1=thr_vic_array(2);
  %m1_enh=m1(2);
  
  %sens_spec_m1_chrN=[m1_enh,vic_m1,chrN];
  sens_spec_m2_chrN=[m2_enh,vic_m2,chrN];% combined between lineages
  %ssm1=sens_spec_m1_chrN;
  ssm2=sens_spec_m2_chrN;
  CAC_opt=max([CAC_opt_ect,CAC_opt_end,CAC_opt_mes]);
end% func  main performance
%===================subf
function [num_enh_1,num_enh_2,num_enh_3]=num_enh_lineage_BDH(BDH1,BDH2,BDH3)

%check for non empty sets of enhs
  if size(BDH1)>0,
    num_enh_1=length(BDH1(:,1));
  else
    num_enh_1=0;
  end
  if size(BDH2)>0,
    num_enh_2=length(BDH2(:,1));
  else
    num_enh_2=0;
  end
  if size(BDH3)>0,
    num_enh_3=length(BDH3(:,1));
  else
    num_enh_3=0;
  end
end% func
%======================subf opt vic
function [ area,ma_ect ] = opt_vic( tp_ect,thr_vic_array )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%tp_ect=tpp(:,1);
  ma_ect=max(tp_ect);
  for i=1:length(tp_ect),% number vic array
      if tp_ect(i)==ma_ect,
          area=thr_vic_array(i);
          break
      end
  end
  %area

end

%======================subf opt vic and CAC
function [ zone,ma_tp_ect,ma_dC_ect,CAC_opt,ind_vic] = opt_vic_CAC( tp_ect,thr_vic_array,CACC_ect,dCC_ect,den_lin_ect)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%tp_ect=tpp(:,1);
%CACC_ect,dCC_ect,den_lin_ect- column vectors length of number vicinities
  ma_dC_ect=max(dCC_ect);
  for i=1:length(tp_ect),% number vic array
      if dCC_ect(i)==ma_dC_ect,
          zone=thr_vic_array(i);% zone of influence
          CAC_opt=CACC_ect(i);
          ind_vic=i;
          ma_tp_ect=tp_ect(i);
          break
      end
  end
  %area

end
%===============================subfunction
function [CAC,tpe,fpe,TPe,FPe,FNe,tpg,fpg,TPg,FPg,FNg]=CAC_simplistic_decisionLineage(prop_LMH,ng_LMH,ne_LMH)

CAC=[];% is 3-vector
%tp=[tp_ect,tp_end,tp_mes];
%fp=[fp_ect,fp_end,fp_mes];
tpe=[];% tpe- enhancers
fpe=[];
tpg=[];% tpg- genes
fpg=[];
%----------------------CAC decision1
        TPe=0;TPg=0;
        FPe=0;FPg=0;
        FNe=0;FNg=0;% CAC is negative but there are enh in high GE thr
 %----------------Ecto
        tpe_ect=0;fpe_ect=0;
        tpg_ect=0;fpg_ect=0;
        CAC_ect=prop_LMH(3,1) - prop_LMH(1,1);
        if CAC_ect > 0,% CAC is positive for H GE, accept Ecto gene-enh as TP
            TPe=TPe+ne_LMH(3,1);
            TPg=TPg+ng_LMH(3,1);
            tpe_ect=ne_LMH(3,1);
            tpg_ect=ng_LMH(3,1);
        else % CAC neg: is positive for L GE, accept Ecto gene-enh as FP: potentially
            FPe=FPe+ne_LMH(1,1);% Low GE
            FPg=FPg+ng_LMH(1,1);
            fpe_ect=ne_LMH(1,1);
            fpg_ect=ng_LMH(1,1);
          
     %-----------------FN: neg CAC and high GE enh detected: doubts??????
            FNe=FNe+ne_LMH(3,1);% high GE
            FNg=FNg+ng_LMH(3,1);    
       end
 
 %--------------------Endo
      tpe_end=0;fpe_end=0;
           tpg_end=0;fpg_end=0;
      CAC_end=prop_LMH(3,2) - prop_LMH(1,2);
      if CAC_end > 0,% CAC is positive, accept Endo gene-enh as TP
            TPe=TPe+ne_LMH(3,2);
            TPg=TPg+ng_LMH(3,2);
            tpe_end=ne_LMH(3,2);
            tpg_end=ng_LMH(3,2);
     
      else % CAC is positive for L GE, accept endo gene-enh as FP
            FPe=FPe+ne_LMH(1,2);% Low GE
            FPg=FPg+ng_LMH(1,2);
            fpe_end=ne_LMH(1,2);
            fpg_end=ng_LMH(1,2);

      %-----------------FN: neg CAC and high GE enh detected
            FNe=FNe+ne_LMH(3,2);% high GE
            FNg=FNg+ng_LMH(3,2);    

     end
    
 %----------------Meso
     tpe_mes=0;fpe_mes=0;
          tpg_mes=0;fpg_mes=0;
     CAC_mes=prop_LMH(3,3) - prop_LMH(1,3);
     if CAC_mes,% CAC is positive, accept Ecto gene-enh as TP
             TPe=TPe+ne_LMH(3,3);
             TPg=TPg+ng_LMH(3,3);
             tpe_mes=ne_LMH(3,3);
             tpg_mes=ng_LMH(3,3);

     else % CAC is positive for L GE, accept Meso gene-enh as FP
             FPe=FPe+ne_LMH(1,3);% Low GE
             FPg=FPg+ng_LMH(1,3);    
             fpe_mes=ne_LMH(1,3);
             fpg_mes=ng_LMH(1,3);
    
     %-----------------FN: neg CAC and high GE enh detected
            FNe=FNe+ne_LMH(3,3);% high GE
            FNg=FNg+ng_LMH(3,3);    

     end
     
     CAC=[CAC_ect,CAC_end,CAC_mes];
     tpe=[tpe_ect,tpe_end,tpe_mes];
     fpe=[fpe_ect,fpe_end,fpe_mes];
     tpg=[tpg_ect,tpg_end,tpg_mes];
     fpg=[fpg_ect,fpg_end,fpg_mes];


end % function
