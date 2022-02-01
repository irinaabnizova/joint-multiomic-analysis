function [ng_LMH,ne_LMH,prop_LMH,TabH1,TabH2,TabH3,cGE1,cGE2,cGE3]=num_enh_lineage_11(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,endN,mesN,indEctN,indEndN,indMesN,vec_EcN,vec_EnN,vec_MN,dar1,dar2,dar3,chrN)
   
 % became -1  becaese of simon_11, because of ind_1 january 22, linea
% already separated by chrN though...

          %-------------------13-00 new
          [TabH1,TabH2,TabM1,TabM2,TabL1,TabL2,BDH1,BDH2,BDM1,BDM2,BDL1,BDL2,cGE1,cGE2,HE_bin_hist,ME_bin_hist,LE_bin_hist,re_rm_12,distH1,distH2,distM1,distM2,distL1,distL2,tar_darH,tar_darM,tar_darL,me12_HML,ncG1,ncG2,ngw1,ngw2,numLMH_12]=simon_plot_distHML_EcEn(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,endN,indEctN,indEndN,vec_EcN,vec_EnN,dar1,dar2,chrN);%,nameH,nameHN,nameL,nameLN);
          numEcEn_LMH=numLMH_12;

          [TabH1,TabH3,TabM1,TabM3,TabL1,TabL3,BDH1,BDH3,BDM1,BDM3,BDL1,BDL3,cGE1,cGE3,HE_bin_hist,ME_bin_hist,LE_bin_hist,re_rm_13,distH1,distH3,distM1,distM3,distL1,distL3,tar_darH,tar_darM,tar_darL,me13_HML,ncG1,ncG3,ngw1,ngw3,numLMH_13]=simon_plot_distHML_EcM(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,mesN,indEctN,indMesN,vec_EcN,vec_MN,dar1,dar3,chrN);%,nameH,nameHN,nameL,nameLN);
          numEcM_LMH=numLMH_13;
          
          %==================13=00 old

         %display('number of genes per high GE threshold and lineage');
         num_gene_EcEnM_LMH=[numEcEn_LMH numEcM_LMH(:,2)];
         num_gene_EcEnM_H=num_gene_EcEnM_LMH(3,:);

         %display('number of Enhancers per high GE threshold and lineage');
         [num_enh_1,num_enh_2,num_enh_3]=num_enh_lineage_BDH(BDH1,BDH2,BDH3);


         num_enh_EcEnM_H=[num_enh_1,num_enh_2,num_enh_3];%[length(BDH1(:,1)),length(BDH2(:,1)),length(BDH3(:,1))]

         [num_enh_1,num_enh_2,num_enh_3]=num_enh_lineage_BDH(BDM1,BDM2,BDM3);
         num_enh_EcEnM_M=[num_enh_1,num_enh_2,num_enh_3];%[length(BDM1(:,1)),length(BDM2(:,1)),length(BDM3(:,1))];
         [num_enh_1,num_enh_2,num_enh_3]=num_enh_lineage_BDH(BDL1,BDL2,BDL3);

         num_enh_EcEnM_L=[num_enh_1,num_enh_2,num_enh_3];%[length(BDL1(:,1)),length(BDL2(:,1)),length(BDL3(:,1))];


%display('number of genes & enhs per three GE thresholds and lineage: L M H');
         num_gene_EcEnM_LMH;
         num_enh_EcEnM_LMH=[num_enh_EcEnM_L;num_enh_EcEnM_M;num_enh_EcEnM_H];

%=================compute proportion enh per gene
        ng_LMH=num_gene_EcEnM_LMH;
        ne_LMH=num_enh_EcEnM_LMH;

        for i=1:3,
            for j=1:3,
                if ng_LMH(i,j)>0,
                   prop_LMH(i,j)=ne_LMH(i,j)/ng_LMH(i,j);
                else
                   prop_LMH(i,j)=0;
                end
            end
        end
        %prop_LMH
end% func
%======================subfunc
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
