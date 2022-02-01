function [TabH1,TabH3,TabM1,TabM3,TabL1,TabL3,BDH1,BDH3,BDM1,BDM3,BDL1,BDL3,cGE1,cGE3,HE_bin_hist,ME_bin_hist,LE_bin_hist,re_rm,distH1,distH3,distM1,distM3,distL1,distL3,tar_darH,tar_darM,tar_darL,me13_HML,ncG1,ncG3,ngw1,ngw3,numHML]=simon_plot_distHML_EcEn(thr_vic,thr10,thrh,t,te,gene_names,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN,dar1,dar3,chrN)

% EcEn is only in the name thrExpr_EcEn_Dif_dist_GEint_ind_1
%%===============hand made histo 1K or 5K: if small enough
%if thr_vic <= thrh,then 1K, esle 5K

%---------------------------------OUTPUTS

%cGE1=[ccE11 GGE1];% count gene Expression. counts divided by nw
%cGE3=[ccE33 GGE3];% count gene Expression

%LE_bin_hist=[bin5;Ncb15;Ncb35]; % 1K-histos of dis2start to TSS- per H M L
%HE_bin_hist=[bin5;Ncb15;Ncb35];
%distH1,distH3- distances bw neighbouring genes: H M L

%tab1,tab3,bbed_dar1,bed_dar3--------------only high Expr,
%re_rm= corr 1 and 3 simlified CAC type (three GE thr)
%me13_HML=[me1_me3_H;me1_me3_H;me1_me3_L];

%-------------------------initiate histo TSS counts
%[TabH1,TabH3,TabM1,TabM3,TabL1,TabL3,BDH1,BDH3,BDM1,BDM3,BDL1,BDL3,cGE1,cGE3,HE_bin_hist,ME_bin_hist,LE_bin_hist,re_rm,distH1,distH3,distM1,distM3,distL1,distL3,tar_darH,tar_darM,tar_darL,me13_HML,ncG1,ncG3,ngw1,ngw3]=simon_plot_1K_distHML_divideNorm_vic(thr_vic,thr10,thrh,t,te,gene_names,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN,dar1,dar3,chrN);

%======================params

thr70=thr_vic;% not band but vicinity here
 
%---------------------initialise
numHML=[];
me13_HML=[];
cGE1=[];cGE3=[];
ncG1=[];ncG3=[];
LE_bin_hist=[];
HE_bin_hist=[];
ME_bin_hist=[];

re_rm=[];
        
%--------------initiate DAR counts
GGE1=[];
ccE1=[];
GGE3=[];
ccE3=[];

ncG1=[];
ncG3=[];
  %==============================high GE bin 3
  for i=3,%:length(t),
    thrExS=t(i);
    thrExE=te(i);%thrExS+step;
    GE_interval=[thrExS,thrExE];

    %d2s =distance2start
     
    [TabH1,TabH3,BDH1,BDH3,tar_darH,d2s1,d2s3,cE1,cE3,GE1,GE3,numH,distH1,distH3]=thrExpr_EcEn_Dif_dist_GEint_ind_2(thr_vic,thrExS,thrExE,gene_names,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN,dar1,dar3,chrN);
    [HE_bin_hist,newh_cGE1,newh_cGE3,Ec_M_H,den_EcM_H,den_vic_EcM_H,nghw1,nghw3,me13_H]=get_features_TSS_scatter(thr_vic,thr10,thr70,thrh,d2s1,d2s3,TabH1,TabH3,distH1,distH3,tar_darH);

  end %i  

     nnumInt_13c(i,:)=numH;%n13;% number genes in the thrExpr interval: Ec and Mes
     TD(:,:,i)=tar_darH;% all genes
     dden_vic_EcM(i,:)=den_vic_EcM_H;
     numEcM(i,:)=Ec_M_H;
     denn_EcM(i,:)=den_EcM_H;
     %-------------------------accumulate
     GGE1=[GGE1;GE1];
     ccE1=[ccE1;cE1];
     GGE3=[GGE3;GE3];
     ccE3=[ccE3;cE3];
     ncG1=[ncG1; newh_cGE1];
     ncG3=[ncG3; newh_cGE3];

     
  for i=2,%:length(t),

    thrExS=t(i);
    thrExE=te(i);%thrExS+step;
    GE_interval=[thrExS,thrExE];
   
      
    [TabM1,TabM3,BDM1,BDM3,tar_darM,d2s1,d2s3,cE1,cE3,GE1,GE3,numM,distM1,distM3]=thrExpr_EcEn_Dif_dist_GEint_ind_2(thr_vic,thrExS,thrExE,gene_names,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN,dar1,dar3,chrN);
    [ME_bin_hist,newm_cGE1,newm_cGE3,Ec_M_M,den_EcM_M,den_vic_EcM_M,ngmw1,ngmw3,me13_M]=get_features_TSS_scatter(thr_vic,thr10,thr70,thrh,d2s1,d2s3,TabM1,TabM3,distM1,distM3,tar_darM);

  end %i  
     nnumInt_13c(i,:)=numM;%n13;% number genes in the thrExpr interval: Ec and Mes
     TD(:,:,i)=tar_darM;% all genes
     dden_vic_EcM(i,:)=den_vic_EcM_M;
     numEcM(i,:)=Ec_M_M;
     denn_EcM(i,:)=den_EcM_M;
     %-------------------------accumulate
     GGE1=[GGE1;GE1];
     ccE1=[ccE1;cE1];
     GGE3=[GGE3;GE3];
     ccE3=[ccE3;cE3];
     ncG1=[ncG1; newh_cGE1];
     ncG3=[ncG3; newh_cGE3];

 
  for i=1,%:length(t),

    thrExS=t(i);
    thrExE=te(i);%thrExS+step;
    GE_interval=[thrExS,thrExE];
      
    [TabL1,TabL3,BDL1,BDL3,tar_darL,d2s1,d2s3,cE1,cE3,GE1,GE3,numL,distL1,distL3]=thrExpr_EcEn_Dif_dist_GEint_ind_2(thr_vic,thrExS,thrExE,gene_names,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN,dar1,dar3,chrN);
    [LE_bin_hist,newl_cGE1,newl_cGE3,Ec_M_L,den_EcM_L,den_vic_EcM_L,nglw1,nglw3,me13_L]=get_features_TSS_scatter(thr_vic,thr10,thr70,thrh,d2s1,d2s3,TabL1,TabL3,distL1,distL3,tar_darL);

  end %i 

     nnumInt_13c(i,:)=numL;%n13;% number genes in the thrExpr interval: Ec and Mes
     TD(:,:,i)=tar_darL;% all genes
     dden_vic_EcM(i,:)=den_vic_EcM_L;
     numEcM(i,:)=Ec_M_L;
     denn_EcM(i,:)=den_EcM_L;
     
     %-------------------------accumulate
     GGE1=[GGE1;GE1];
     ccE1=[ccE1;cE1];
     GGE3=[GGE3;GE3];
     ccE3=[ccE3;cE3];
     ncG1=[ncG1; newh_cGE1];
     ncG3=[ncG3; newh_cGE3];

%=================================================
  % me1_me3_H=[mean(distH1),mean(distH3)];
   
   %[newh_cGE1,newh_cGE3,nghw1,nghw3]=normalise_countEnh(TabH1,TabH3);
    %bed_dar1 =BD
    %       1.00   14880497.00   14880597.00         36.00   14918862.00   36683183.00             0          1.08
    %      1.00   14920097.00   14920197.00         36.00   14918862.00   36683183.00             0          1.08
    %      1.00   14947897.00   14947997.00         36.00   14918862.00   36683183.00             0          1.08
    %      1.00   14962297.00   14962397.00         36.00   14918862.00   36683183.00             0          1.08
     
             
    %-----------------output for HE genes
    %--------------------'interaction' matrix (other way around)
    %bbed_dar1=BDH1 = bed_dar1;
    %bbed_dar3=bed_dar3;
       
       %-------------------------ACCUMULATE i (H M L GE data)
        
       %   [newh_cGE1,newh_cGE3]=normalise_countEnh(TabH1,TabH3);
 
       %ncG1=[ncG1; newh_cGE1];
       %ncG3=[ncG3; newh_cGE3];
     
  
       %-------------------------accumulated across H M L


       ngw1=nghw1+ngmw1+nglw1;
       ngw3=nghw3+ngmw3+nglw3;
       
       %=====================CAC out of three GE

       [rEc_d pEc]=corrcoef(t,denn_EcM(:,1));
       %pe=pEc(1,2);
       [rM_d pM]=corrcoef(t,denn_EcM(:,2));
       %pm=pM(1,2);
       %[rEn_d pEn]=corrcoef(t,denn_EcEn(:,2));
       %-----------------------correlations for this vicinity: not enough of
        %points, only three binGE intervals (kind of CAC simple1!
      rEcv=rEc_d(1,2);
      %rEnv=rEn_d(1,2);
      rMv=rM_d(1,2);
      re_rm=[rEcv,rMv];

       % ----------------------more compact outputs
      cGE1=[ccE1 GGE1];% count gene Expression
      cGE3=[ccE3 GGE3];% count gene Expression

      me13_HML=[me13_H;me13_M;me13_L];% mean distance bw DEGs- we know differently
      numHML=nnumInt_13c;

end 
%============================subfunctions
function [Tab_G1_DAR1,Tab_G3_DAR3,bed_dar1,bed_dar3,tar_dar,d2s1,d2s3,cE1,cE3,GE1,GE3,num,dist1,dist3]=thrExpr_EcEn_Dif_dist_GEint_ind_2(thr_vic,thrExS,thrExE,gene_names,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN,dar1,dar3,chrN)


%----------------will output BDH dar/enh index and gene GE, 11 July 2020


%INPUTS: sh be matched dar_cm and dar_mc
%1. =====================Diff and then  per chromosome

%=============================OUTPUTS
%nCEcM=num_Comm_dEc_dM;- all chroms Diff
%numInt_13, numInt_13c=num1_3;
%sl=[sl1(1),sl3(1)];
    %bed_chat_dar3= [G3_DAR3(:,1)   DAR3_genes  G3_DAR3(:,5)] ;      
    %bed_chat_dar1= [G1_DAR1(:,1)   DAR1_genes  G1_DAR1(:,5)] ;    
    %bcd3=bed_chat_dar3;
    %bcd1=bed_chat_dar1;
  %-----------initiate empty outputs

Tab_G1_DAR1=[];
Tab_G3_DAR3=[];
bed_dar1=[];
bed_dar3=[];

tar_dar=[[0,0]' [0,0]'];
%sl=[];,d2s1,d2s3,cE1,cE3,GE1,GE3,num,
bcd1=[];
bcd3=[];

dist1=[];dist3=[];

sl=[];
d2s1=[];d2s3=[];
cE1=[];cE3=[];GE1=[];GE3=[];
num=[0,0];
tdo=[[1,0];[1,0]];

%display(' this vicinity and ExprThr: linked genes and DARs');
%int_Exp=[thrExS,thrExE]


%1.2======================Genes per chromosome chrN
[num1_3,G1In_exp,g1In_ind,G3In_exp,g3In_ind]=DifExpr_interval(thrExS,thrExE,geneEc_chrN,geneM_chrN,ind_geneEcN,ind_geneMN,vec_EcN,vec_MN);
num=num1_3;

%==========================per lineage 

  if min(size(G1In_exp))>0,
    
     [Tab_G1_DAR1,G1_DAR1,names_linked_G1,DAR1_genes,nonRedun1,tar_dar_1,d2s1,GE1,countEnh1,sl1,pos_dir1]=link_EnhInd_tarGene_any_1(thr_vic,G1In_exp,g1In_ind,G1In_exp(:,5),G1In_exp(:,6),gene_names,dar1,chrN);
 
     %[gene_names_DAR1,Tab_G1_DAR1,G1_DAR1,names_linked_G1,DAR1_genes,DAR1_genes_ou,TGD_in1,TGD_ou1,nonRedun1,tar_dar_1,d2s1,GE1,countEnh1,sl1,dio1,pos_dir1]=link_Enh_tarGene_dist_en(thr_vic,G1In_exp,g1In_ind,G1In_exp(:,5),G1In_exp(:,7),gene_names,dar1,chrN);
     %title(['GE1 vs DAR1 count, slope=',num2str(sl1(1))]);
  
    [dist1]=distances_analysis_np(pos_dir1);
    me_dist1=mean(dist1);
    cE1=countEnh1;

  else
    display([' no data for GE1 this thrE=',num2str(thrExS)]);
  end

%---------------------

%display(' Mes for this thrEx');
  if min(size(G3In_exp))>0,
   
      [Tab_G3_DAR3,G3_DAR3,names_linked_G3,DAR3_genes,nonRedun3,tar_dar_3,d2s3,GE3,countEnh3,sl3,pos_dir3]=link_EnhInd_tarGene_any_1(thr_vic,G3In_exp,g3In_ind,G3In_exp(:,6),G3In_exp(:,5),gene_names,dar3,chrN);

      %[gene_names_DAR3,Tab_G3_DAR3,G3_DAR3,names_linked_G3,DAR3_genes,DAR3_genes_ou,TGD_in3,TGD_ou3,nonRedun3,tar_dar_3,d2s3,GE3,countEnh3,sl3,dio3,pos_dir3]=link_Enh_tarGene_dist_en(thr_vic,G3In_exp,g3In_ind,G3In_exp(:,7),G3In_exp(:,5),gene_names,dar3,chrN);
      %title(['GE3 vs DAR3 count, slope=',num2str(sl3(1))]);
      [dist3]=distances_analysis_np(pos_dir3);
      me_dist3=mean(dist3);
      cE3=countEnh3;

  else
    
    display([' no data for GE2 this thrE=',num2str(thrExS)]);
  end

   %-------------------if both Tabs are non-empty

  if size(Tab_G1_DAR1)>0 & size(Tab_G3_DAR3) >0,

   %tar_dar11=[length(gene_names_DAR1),length(names_linked_G1)];
   %tar_dar33=[length(gene_names_DAR3),length(names_linked_G3)];

   %sl=[sl1,sl3];
   cE1=countEnh1;
   cE3=countEnh3;

   tar_dar=[tar_dar_1;tar_dar_3];
   %si=size(TGD_ou1)
   %if size(TGD_ou1)>0 & size(TGD_ou3)>0,
   %tar_dar_out=[[height(TGD_ou1),length(DAR1_genes_ou(:,1))];[height(TGD_ou3),length(DAR3_genes_ou(:,1))]];
   %tdo= tar_dar_out;
   %end

   
   %===============Final bed dar =BD file: 9 fields from 12 July (DARi_gene
   %las col is dar index
   
   %BD= chr dar_st dar_end gene_ind TSS TES dir  GE'); 8 fields

    bed_dar3= [G3_DAR3(:,1)   DAR3_genes  G3_DAR3(:,2:5)] ;      
    bed_dar1= [G1_DAR1(:,1)   DAR1_genes  G1_DAR1(:,2:5)] ; 
    
    
    %-------------------intermediate files
    
    %G3_DAR3 =

        %  3.00   14863538.00   14872351.00          1.00          4.95          1.10       2450.00          2.00
        %  3.00   53041528.00   53261679.00          1.00          5.95          2.33       2564.00          1.00
        %  3.00   53041528.00   53261679.00          1.00          5.95          2.33       2564.00          1.00
        %  3.00   58692400.00   58692400.00             0          0.34         -2.12       2597.00          1.00
        %  3.00   58692400.00   58692400.00             0          0.34         -2.12       2597.00          1.00
      

    %bed_chat_dar3= [G3_DAR3(:,1)   DAR3_genes  G3_DAR3(:,5)] ;      
    %bed_chat_dar1= [G1_DAR1(:,1)   DAR1_genes  G1_DAR1(:,5)] ;    
    %bcd3=bed_chat_dar3;
    %bcd1=bed_chat_dar1;
    
  else
    display(['no data GE1 or GE2 for these thrE and vicinity=',num2str(thr_vic)]);
  end
end % function



%=========================subfunction 2
function [HE_bin_hist,newh_cGE1,newh_cGE3,Ec_M,den_EcM,den_vic_EcM,nghw1,nghw3,me1_me3_H]=get_features_TSS_scatter(thr_vic,thr10,thr70,thrh,d2s1,d2s3,TabH1,TabH3,distH1,distH3,tar_darH)

%-------------------outputs
%HE_bin_hist (HE_bin_hist=[bin5;Ncb15w;Ncb35w];- three horison vectors
%newh_cGE1,newh_cGE3,
%nghw1,nghw3
%Ec_M,den_EcM,den_vic_EcM,

    me1_me3_H=[mean(distH1),mean(distH3)];
   
   [newh_cGE1,newh_cGE3,nghw1,nghw3]=normalise_countEnh(TabH1,TabH3);
    %bed_dar1 =BD
    %       1.00   14880497.00   14880597.00         36.00   14918862.00   36683183.00             0          1.08
    %      1.00   14920097.00   14920197.00         36.00   14918862.00   36683183.00             0          1.08
    %      1.00   14947897.00   14947997.00         36.00   14918862.00   36683183.00             0          1.08
    %      1.00   14962297.00   14962397.00         36.00   14918862.00   36683183.00             0          1.08
   
     

   Ec_M=[tar_darH(1,1),tar_darH(2,1)];

   [den_EcM]=density_interval(thr10,thr70,d2s1,d2s3,Ec_M);%  vic or band!

   den_vic_EcM=[tar_darH(1,2)/tar_darH(1,1),tar_darH(2,2)/tar_darH(2,1)];

   %===============hand made histo 1K or 5K
   if thr_vic <= thrh,

    %-----------------------from May2020
   [bin5, count_bin15,Ncb15,Ncb15w]=histo1K_NormVic(thr_vic,d2s1,Ec_M(1),nghw1);
   [bin5, count_bin35,Ncb35,Ncb35w]=histo1K_NormVic(thr_vic,d2s3,Ec_M(2),nghw3);

   else
    
   [bin5, count_bin15,Ncb15,Ncb15w]=histo5K_normVic(thr_vic,d2s1,Ec_M(1),nghw1);
   [bin5, count_bin35,Ncb35,Ncb35w]=histo5K_normVic(thr_vic,d2s3,Ec_M(2),nghw3);
  
   end

   HE_bin_hist=[bin5;Ncb15w;Ncb35w];

   
end % function
%==========================subfunc3
%function [newh_cGE1,newh_cGE3,nghw1,nghw3]=normalise_countEnh(TabH1,TabH3);
function [new_cGE1,new_cGE3,ngw1,ngw3]=normalise_countEnh(TabH1,TabH3)

% divides countEnh by number of genes in vicinity of main gene
%Input
%TABH1(1:30,:) here SE genes
%                   index                      pos_dir                     countEnh    GE1     GE2      NW 
 %                    ______    _________________________________________    ________    ____    ____    ____
%
%    Tcea1              6.00     4857814.00     4897909.00           1.00     3.00       5.31    5.19    3.00
%    Rrs1              16.00     9545408.00     9547454.00           1.00     2.00       5.49    5.46    1.00
%    Cops5             23.00    10038127.00    10024602.00           0.00     3.00       5.84    5.89    1.00
%    Tram1             31.00    13589864.00    13564702.00           0.00     6.00       5.67    5.42    1.00
   


%------------outputs
%new_cGE1=[new_cenh1' GE1];
 %new_cGE3=[new_cenh3' GE3];
 
 %----------------initialis
 ngw1=0;% summ num geens in the vicinity
 ngw3=0;
 
 
 new_cGE1=[];
 new_cGE3=[];

if size(TabH1)>0,
    
nw1=TabH1{:,{'NW'}};
 nenh1=TabH1{:,{'countEnh'}};
 for i=1:length(nw1),
      new_cenh1(i)= nenh1(i)/nw1(i);
  end
GE1=TabH1{:,{'GE1'}};
new_cGE1=[new_cenh1' GE1];

ngw1=sum(nw1);
end



%=====================Mes
 if size(TabH3)>0,
     
  nw3=TabH3{:,{'NW'}};
  nenh3=TabH3{:,{'countEnh'}};
  
   for i=1:length(nw3),
     new_cenh3(i)= nenh3(i)/nw3(i);
   end
  
   GE3=TabH3{:,{'GE1'}};
   new_cGE3=[new_cenh3' GE3];
   
   ngw3=sum(nw3);
 end

end %function
%=========================f4
%[den_EcM]=density_interval(thr10,thr70,d2s1,d2s3,Ec_M);%  vic or band!
function [den_EcM]=density_interval(thr10,thr70,dis2startE,dis2startM,Ec_M)

 %upper threshold is max of distances computed
 thrMax=max(max(dis2startE),max(dis2startM));
 a=thrMax-thr70;
 ne=0;
 for i=1:length(dis2startE),
     if abs(dis2startE(i))>thr10,
         ne=ne+1;
     end
 end
 denEc=ne/Ec_M(1);
 
 nm=0;
 for i=1:length(dis2startM),
     if abs(dis2startM(i))>thr10,
         nm=nm+1;
     end
 end
 denM=nm/Ec_M(2);
 
 %thrExS
 den_EcM=[denEc,denM];

end % function
%=======================f5
%[bin5, count_bin15,Ncb15,Ncb15w]=histo1K_NormVic(thr_vic,d2s1,Ec_M(1),nghw1);
function [bin, count_bin,Ncb,Ncbv]=histo1K_NormVic(thr_vic,d2s1,num_genes,ngw)

 %compute  counts (count_bin) of DAR starts in 1KB up/downstream TSSs from d3s 
 %Ncb= Normalised count bin

 % ngw total number of genes, includingthose in the vicinity of TSS: it
 % it might be >num_genes
 
 %Ncbv=count_bin/ngw; Normalised count bin vicinity
 
 bin=-thr_vic:1000:thr_vic;
 
 %------initialiise per bin counts
 
   for j=1:length(bin),
     count_bin(j)=0;
   end
 
   sd2s=sort(d2s1);
 
   i=1;

  for j=1:length(count_bin)-1,
    
    while (i<=length(sd2s) & sd2s(i) <= bin(j+1) ),
        count_bin(j)= count_bin(j)+ 1;
         i=i+1;
    end
  end
 %count_bin to normalise per gene

 Ncb=count_bin/num_genes;

Ncbv=count_bin/ngw;

end % function
%=========================subfun6
%[bin5, count_bin15,Ncb15,Ncb15w]=histo5K_normVic(thr_vic,d2s1,Ec_M(1),nghw1);
function [bin, count_bin,Ncb,Ncbv]=histo5K_normVic(thr_vic,d2s1,num_genes,ngw)

 %compute  counts (count_bin) of DAR starts in 1KB up/downstream TSSs from d3s 
 %Ncb= Normalised count bin

 % ngw total number of genes, includingthose in the vicinity of TSS: it
 % it might be >num_genes
 
 %Ncbv=count_bin/ngw; Normalised count bin vicinity
 
 bin=-thr_vic:5000:thr_vic;
 
 %------initialiise per bin counts
 
 for j=1:length(bin),
     count_bin(j)=0;
 end
 
 sd2s=sort(d2s1);
 
 i=1;

for j=1:length(count_bin)-1,
    
    while (i<=length(sd2s) & sd2s(i) <= bin(j+1) ),
        count_bin(j)= count_bin(j)+ 1;
         i=i+1;
    end
end
%count_bin to normalise per gene

Ncb=count_bin/num_genes;

Ncbv=count_bin/ngw;

end %function

%===============func7
function [numInt_EcM,geneEcIn_exp,EcIn_ind,geneMIn_exp,MIn_ind]=DifExpr_interval(thrExS,thrExE,DEcgene,DMgene,indEc,indM,vecEc,vecM)
%separates genes within Expr interval, 
%computes their indexes and number within interval: for ANY two layers  

%=====================INPUTS

%DEcgene,DMgene,indEc,indM,vecEc,vecM)= might be for ANY (here written Ec mes)
%  two layers, but 

% main file, GE and ind should match corresponding layer

%=======================OUTPUTS
numInt_EcM=[0,0];%number of genes in a given interval (thrExS,thrExE)


geneEcIn_exp=[];
EcIn_ind=[];
geneMIn_exp=[];
MIn_ind=[];

% Ec vs mes (7.5) or EPI vs TB if 6.5
%L13= Lineages 1 vs 3 (Ect Mes or EPI TB)

%Input: index Dif genes

%indEc, indM is per chrom gene index

%-------------------OUTput indexes for Intevals Expression genes

%Chromosome	Start	End	Probe.Strand	E6.5_EPI	E6.5_PS	E6.5_TB	E6.5_VE
%1	3205901	3671498	0  -10.13058	0.6720925	-7.3292418	-4.0362945
%1	4343507	4360314	0  -10.43162	-6.6878405	-8.782284	-8.074299
%1	4490928	4496413	0  -2.7995515	-5.4365945	-7.531038	4.5727186

%---------------------gets interval for Expression level: No PlotsOutputs
%and no SN!!!!


%-------------------separate high low Expressed
%Chromosome	Start	End	Probe.Strand	E7.5_Ectoderm	E7.5_Endoderm	E7.5_Mesoderm
%1	3205901	3671498 0 -4.35045	-2.3982441	-3.20762

%thrExEx=3;%was -3

%geneM=gene(:,7);
%geneEn=gene(:,6);
%geneEc=gene(:,5);


    
 %vecEc= DEcgene(:,5);
%1===================================Ect DiffEcM_E=DEcgene
%display('Diff Ect genes within expr interval'); 
[geneEcIn_exp,geneEcH_exp,geneEcL_exp,EcIn_ind]=separate_expression_interval_ind(thrExS, thrExE,DEcgene,vecEc,indEc);
if min(size(geneEcIn_exp))>0,   
numInEc=length(geneEcIn_exp(:,1));
       %numInt_EcM=[numInEc,numInM];
else
 numInEc=0;
end
    %numInEc
 %vecM=DMgene(:,7);% vec3!!!!!

%2===========================Meso = DiffEcM_M=DMgene
%display('Diff Mes genes within expr interval'); 
[geneMIn_exp,geneMH_exp,geneML_exp,MIn_ind]=separate_expression_interval_ind(thrExS, thrExE,DMgene,vecM,indM);
if min(size(geneMIn_exp))>0,
 numInM=length(geneMIn_exp(:,1));
else
  numInM=0;
end
%numInM


numInt_EcM=[numInEc,numInM];
end% func
%=================func8
function [geneMIn_exp,geneMH_exp,geneML_exp,geneMIn_ind]=separate_expression_interval_ind(thrS,thrE, gene_expr,geneM,gene_ind)
%sepatates gene expression values within given thr interaval
%gene_expr
%Chr Start	End	   Strand	E7.5_Ectoderm	E7.5_Endoderm	E7.5_Mesoderm
%1	3205901	3671498	0	04.35045	02.3982441	03.20762
%1	4343507	4360314	0	09.973417	08.495143	03.2783616

%geneM=gene_expr(:,7);

%-----in case of all croms, gene ind sh be ind=1:length(gene)

geneMH_exp=[];
geneML_exp=[];
geneMIn_exp=[];
geneMIn_ind=[];

thr_start=thrS;
thr_end=thrE;

ci=0;
ch=0;clo=0;
  for i=1:length(geneM),
    
    if geneM(i)<= thrE && geneM(i)>=thrS,
    ci=ci+1;
    geneMIn_exp(ci,:)=gene_expr(i,:);
    geneMIn_ind(ci)=gene_ind(i);
    end
    
    if geneM(i)>thrE, 
    ch=ch+1;
    geneMH_exp(ch,:)=gene_expr(i,:);
    end
    
    if geneM(i)<thrS, 
    clo=clo+1;
    geneML_exp(clo,:)=gene_expr(i,:);
    end

  end
%
  if ci < 1,
    display(' no genes within this interval');
  end
end% fun


%====================subf9
function [Tab_G1_DAR1,genes_DAR1,names_linked_genes,DAR1_genes,nonRedun1,num_NZtar_Enh,dis2start,GE1,countEnh,sl,pos_dirP]=link_EnhInd_tarGene_any_1(thr_vic,targett,ind_target,exprr1,exprr2,gene_names,Enh,chrN)
   %----------------------any means: just compute n-w in cluster fat, not
   %discard based on tobe the only1 in the vicinity
   %Links diff expressd genes (expr1 vs expr2) to their diff accessible enhancers/DARs

    %dense_in_out
   %=======================e.g. Ecto vs Meso (1 vs 2)

   num_NZtar_Enh=[];
   Tab_G1_DAR1=[];

 
  genes_DAR1=[];
  names_linked_genes={};
  DAR1_genes=[];
  nonRedun1=[];
  countEnh=[];
  TargetEnh=[];
  EnhTarget=[];
  count_Enh_target=[];
  dis2start=[];% distances for enhancers relative to nearest gene start
  indGene_targetEnh=[];
 
  perc=[];
  posEnh=[];
  len=[];
  pos_tar_enh=[];
  
  %------------------------indexes
  ind_Enh=[];
  
  ch=[];
  new_gene=[];
  
  %------------------------------possible distances
  pos_dirP=[];

 %------------------------- if non-empty data

 if min(size(targett))>0,

    pos_dir=targett(:,2:4);

  %===============filter out too close genes: chose one with bigger GE
%for a current lineage!
  %======================cluster near-by genes, and reduce their amount if too
%close ones

 [intervals,n_g,new_gene,clus_st,clusEc,clusM,clusInd,clusDir,clus_en]=cluster_fat_pos_difExp_en(pos_dir,exprr1,exprr2,ind_target,thr_vic);

 % here the max Expression in 'new_gene'  for the first input GE vector, ,exprr1
  % new_gene=[start_position' en' direction' gEx1' gEx2' GInd' n_g'];

%===============================new_gene filtered values

  %let this column be chrN  
  for i=1:length(n_g),
    ch(i)=chrN;
  end

  target=[ch' new_gene];% just extra first column - for further computing of positions.
  expr1=new_gene(:,4);
  expr2=new_gene(:,5);
  ind_target=new_gene(:,6);
  
%==========================vicinity around the starts of targets

%-----------------------------if non-empty target data
  
  if max(size(target))>0,
    
    start_target=target(:,2);
    posEnh=Enh(:,1);%  ends are colunm 2
    end_target=target(:,3);
    dir=target(:,4);
    %ind_target=target(:,6);
      
    %========================--initiate number of enhancers around each target start position
    for j=1:length(start_target),
        count_Enh_target(j)=0;
    end
    
    %==========================count all enhancers around all targets
    cenh=0;% count all enhancers around all targets
   ni=0;no=0;
    for j=1:length(start_target),% all given genes
       for i=1:length(posEnh), % all given dars TSS
           
           if abs(posEnh(i)-start_target(j)) < thr_vic,
               
            count_Enh_target(j)=count_Enh_target(j)+1;% per target: count enhancers located in its thr-neighbourhood
            
            cenh=cenh+1;% count all dars/enhancers located in the thr-neibourhood of chrN target targets
            
            %----------------features of linked enhancers and targets
            dis2start(cenh)=posEnh(i)-start_target(j);
            
            %======================Enh: for Main2 output1  %Enh=DARs here within targets vicinity
            EnhTarget(cenh,:)=Enh(i,:);% set of Enhancer around target
            ind_Enh(cenh)=i;% index in the Enh (do I need it?) Y 12 July 2020
            
            %======================targets; %=============for Main2 output2
            TargetEnh(cenh,:)=target(j,:);
            indGene_targetEnh(cenh)=ind_target(j);% Enh-related indexes in the file of targets
            
            
            len(cenh)=abs(end_target(j)-start_target(j));% for inside normalization
            direc(cenh)=dir(j);
            flank(cenh)=len(cenh)-10;% flank one length of a gene: two lengths?...
            
            exprT(cenh)=expr1(j);% strong dif Expr 
            exprTM(cenh)=expr2(j);% meso expr for target genes
            
          
            %------diff?
           end % if
       end % i
   end % j

   numEnh_involved=cenh;
     
   %=====================select Non-redundant genes with at least one enh around
   %from count_Enh_target
            
      [count_Enh_targetP,indP,exprP1,exprP2,pos_dirP, n_w_gene,lenP]=select_nonRedun(count_Enh_target,expr1,expr2,ind_target,target);
      num_nonZ_target=length(lenP);

      %display(' number of targets and Enhancers involved');% tar dar?
      num_NZtar_Enh=[num_nonZ_target,numEnh_involved];
      names_linked_genes=gene_names(indGene_targetEnh);

    %=====================MAKE NON redundant Table
    %Tab_G1_DAR1=table(index,pos_dir,countEnh,GE1,GE2,NW,'RowNames',gene_names_DAR1);
    %nonRedun1=[index pos_dir countEnh GE1 GE2 NW];
    
    [Tab_G1_DAR1,nonRedun1,GE1,countEnh,sl]=make_table_gene_darCountNP(count_Enh_targetP,indP,exprP1,exprP2,pos_dirP, n_w_gene,gene_names);
    %GE1-vector of gene expressions Ect
    %[Tab_G1_DAR1_chrN,nonRedun1,countEnh,sl,r]=make_table_gene_enhCount(count_Enh_targetP,indP,exprP1,exprP2,exprP3,pos_dirP, n_w_gene,gene_names);

  genes_DAR1=[TargetEnh];
  geneDAR=genes_DAR1;

  %======================Main2 output
   if size(EnhTarget)>0,
    DAR1_genes=[EnhTarget(:,1:2) ind_Enh' indGene_targetEnh'];
   end
    
   else
    display('-----------no DARs/enhs around these targets');
  end % c>0
 
 else
     display('-----------empty targett');

 end % if max(size(target))>0,%    main1 if non empty data target
end

%======suf 10
function [intervals,n_g,new_gene,clus_st,clusEc,clusM,clusInd,clusDir,clus_en]=cluster_fat_pos_difExp_en(pos_dir,exprEc,exprM,ind,thr_vic);

%----------------with proper ends of rev genes!

% here the max Expression in 'new_gene'  for the first input GE vector, ,expr
% new_gene=[start_position' en' direction' gEc' gM' GInd' n_g'];

%cluster is represented by one gene with highest GE of exprEct, gEc (first input
%vector


% Clusters genes due to the closenest of their ordered positions (done by
% chromosome!). Depth of clustering is by thr_vicinity

% computes number of genes in a cluster, n_g

%  one gene, giving max to the first Expr input, is selected to represent a cluster 

  %------------------------Inputs: three column matrix
%pos_dir(1:4,:) genes: may be filtered genes
%    pos_start   pos-end         dir
 %   3678115.00    3813122.00          1.00
%    3835665.00    3831334.00             0
 %   3872105.00    3870657.00             0

 %-----------------ind: 
 %   a.index of original gene list, e.g. at elas chromosome, diff etc
 %   b. just ind=1:length(pos_dir(:,1))- not here!
 
 %-------------exprEc=GE (per each gene), vector, exprM-expr vec Mes
 
 
%---------OUTputs

%new_gene=[start_position' en' direction' gEc' gM' GInd' n_g'];

%---------intervals: start- end of gene or gene cluster (if near by)
%    n_g= # genes in a cluster
  
%where the first column is pos-vector of sorted (increasing) genome positions

%thr_vic=depth of clustering-------------     
%------------------output

intervals=[];
n_g=[];
clus_st=[];
clus_en=[];
%difEcM=exprEc-exprM;%  need any?....
pos_st=pos_dir(:,1);
pos_en=pos_dir(:,2);
dir=pos_dir(:,3);
%%%%%%%%cluster growth until next position is less than thr_vic 
    
    
     k1=1; % # of intervals in a cluster
     k=1;% # of cluster
     
     clus_st(k,k1)=pos_st(1);
     clus_en(k,k1)=pos_en(1);
     
     % what are the features in the curr cluster: Expr level Ect mes, ind,
     % dir
     
     clusEc(k,k1)=exprEc(1);
     clusM(k,k1)=exprM(1);
     clusInd(k,k1)=ind(1);
     clusDir(k,k1)=dir(1);
     
     %0. ------------------ initialise the first values
     
     n_g(k)=1; % # of intervals in k-th cluster
     
     st(k)=pos_st(1);% start of first cluster 
     ppos(k)=pos_st(1);% positions skipping clusters, almost starts of clusters
   
    
    %1 ======================cluster by start positions
     
   if length(pos_st) > 1,% list pos more than one value
       
      for i=1:length(pos_st)-1;% 4 March 19 ; 29 Dec 2017 change
        if abs(pos_st(i+1)-pos_st(i))<=thr_vic,% add pos_st to kth cluster
        k1=k1+1;
        clus_st(k,k1)=pos_st(i+1);
        clus_en(k,k1)=pos_en(i+1);
        clusEc(k,k1)=exprEc(i+1);
        clusM(k,k1)=exprM(i+1);
        clusInd(k,k1)=ind(i+1);
        clusDir(k,k1)=dir(i+1);
  
        n_g(k)=k1;
        
        %en(k)=pos_en(i);%pos_st(i+1)+thr_vic; % end of k-th cluster
        else
        %en(k)=max(pos_en(i),pos_en(i-1));%+ending; % end of k-th cluster:CHANGE!
        %en(k)=max(pos_st(i)+5000,pos_en(i));%+ending; % end of k-th cluster:CHANGE!
        en(k)=pos_en(i);%max(pos_en(i),pos_st(i));% end of last gene in a cluster; if rev, then st(i)
        
        %1.1 =======================start new cluster
        ppos(k)=pos_st(i);
        k=k+1;
        k1=1;
        
        clus_st(k,k1)=pos_st(i+1);
        clus_en(k,k1)=pos_en(i+1);
        clusEc(k,k1)=exprEc(i+1);
         clusM(k,k1)=exprM(i+1);
        clusInd(k,k1)=ind(i+1);
        clusDir(k,k1)=dir(i+1);
        st(k)=pos_st(i+1);
        n_g(k)=k1;
        %n_e(k)=num_ev(i+1);       
        end % current cluster formation
     end% for i
    
       % the last cluster end!
      % en(k)=pos_en(length(pos_st));
     en(k)=max(pos_en(length(pos_st)),pos_st(length(pos_st)));%+ending;%  if last gene isrev
    
     ppos(length(st))=pos_st(length(pos_st));
  
  else
         en(1)=pos_en(1);% one item in a list of positions (????)
  end % if length(pos)>1
     
  %------------------------------CLUSTERS are:
   intervals=[st' en'];
   
   
   %2=================get most Ect expressed gene in a cluster...(???? is it
   %OK?)

   %----------------------pick up one gene per cluster as with max GE for Ec
   %here
    gEc=[];GInd=[];start_position=[];direction=[];
  for j=1:length(n_g),
      %------------compute number non-zero Expressions: 
      %-------find max of it, its ind and pos
      exppEc=clusEc(j,:);
      exppM=clusM(j,:);
      indd=clusInd(j,:);
      poss=clus_st(j,:);
      poss_en=clus_en(j,:);
      dirr=clusDir(j,:);
      %[maEx,index,start_pos,direct]=maxEx_non_zero_elements(expp,indd,poss,dirr);
      % what to do if all are zero elements?
      [maEx,expM,index,start_pos,end_pos,direct]=maxEx_non_zero_elements_dif_en(exppEc,exppM,indd,poss,poss_en,dirr);

      gEc(j)=maEx;
      gM(j)=expM;
      GInd(j)=index;
      start_position(j)=start_pos;
      end_position(j)=end_pos;
      direction(j)=direct;
    
  end  
    
%[intervals(1:10,:) n_g(1:10)' start_position' gE'  GInd' ]    
   
   new_gene=[start_position' end_position' direction' gEc' gM' GInd' n_g'];

end %func
%=================subf11
function [count_Enh_targetP,indP,exprP1,exprP2,pos_dirP, n_w_gene,lenP]=select_nonRedun(count_Enh_target,expr1,expr2,ind_target,target)

%=====================select Non-redundant genes with at least one enh around
   %from count_Enh_target
count_Enh_targetP=[];
indP=[];exprP1=[];
exprP2=[];pos_dirP=[]; n_w_gene=[];lenP=[];
c=0;
  for j=1:length(count_Enh_target),
    if count_Enh_target(j)>0,
        c=c+1;
        count_Enh_targetP(c)=count_Enh_target(j);
        exprP1(c)=expr1(j);
        exprP2(c)=expr2(j);
        indP(c)=ind_target(j);%  !!! important
        %difEcM_P(c)=difEcM(j);
        pos_dirP(c,:)=target(j,2:4);
        n_w_gene(c)=target(j,8);%  to check"=number genes around
        lenP(c)=pos_dirP(c,2)-pos_dirP(c,1);
    end
  end
num_nonZ_target=c;
% numEnh_involved=cenh;
  if c==0,
    display(' no linked henes enh this band and thrEx');
  end
end % func

%===================subf12
function [maEx,exppM,index,start_pos,end_pos,direction]=maxEx_non_zero_elements_dif_en(expEc,expM,ind,pos,pos_en,dir)

%given three vectors with concordant number of zeros, find those corr
%elements | giving max to Expr Ecto, not meso
    %expEc=clusEx(j,:);
   % ind=clusInd(j,:);
    %pos=clus(j,:);
    %dir =clusDir(j,:)
    
    maEx=0;
    exppM=0;%
    index=0;
    start_pos=0;
    end_pos=0;
    direction=0;
    
    if size(expEc)>0 & expEc(:,1)~=0,% non empty and non-zero
%1--------find nonZ expEc
  nz=0;expNZ=[];
  for i=1:length(expEc),
    if expEc(i)~= 0,
        nz=nz+1;
        expNZ(nz)=expEc(i);
    end
  end

 maEx=max(expNZ); 
 
 %find what element is it
  for i=1:length(expEc),
     if expEc(i)==maEx,
         jma=i;
     end
  end
 
 index=ind(jma);
 start_pos=pos(jma);
 direction=dir(jma);
 exppM=expM(jma);
 end_pos=pos_en(jma);%------------------16 Feb!!!
 
    %end
    
    else
        index=0;
    end
    %index
 
end

%==================subf 13
function [tab_gene_DAR1_chrN,nonRedun1,GE1,countEnh,sl]=make_table_gene_darCountNP(count_Enh_targetP,indP,exprP1,exprP2,pos_dirP, n_w_gene,gene_names)

% --make table of genes linked to DARe: 
% gene_name gene_ind posSt posEnd dir Expr difEcM_P #Dars 

%inputs: results of [count_Enh_targetP,indP,exprP1,exprP2,pos_dirP, n_w_gene,lenP]=select_nonRedun(count_Enh_target,expr1,expr2,ind_target,target);

nonRedun1=[];
tab_gene_DAR1_chrN=[];

sl=0;
GE1=[];
countEnh=[];

%==========================compute their indexes and NAMES

  if length(indP)>0,
    
    %display('-----------there are DARs around these targets!');
     gene_names_DAR1=gene_names(indP);

     %TABLE ===========make Non-redundand gene table with  counts of DARs
  
    index=indP';
    pos_dir=pos_dirP;
    countEnh=count_Enh_targetP';
    GE1=exprP1';
    GE2=exprP2';
    NW=n_w_gene';
  %DGE=difEcM_P';
  %tab_gene_DAR1_60K_chr4=table(index,pos_dir,countEnh,GE,DGE,'RowNames',gene_names_DAR1);

   %=====================NON redundant Table
    tab_gene_DAR1_chrN=table(index,pos_dir,countEnh,GE1,GE2,NW,'RowNames',gene_names_DAR1);

    nonRedun1=[index pos_dir countEnh GE1 GE2 NW];
  
  
%figure;
%subplot(2,1,1);
%plot(GE1,countEnh,'pr');
%grid;%title('GE vs DAR count');
%xlabel('GE of layer');ylabel('#DARs');
%hold on;
  [inter,slope,ydata22]=slope_int2(GE1,countEnh);
%plot(GE1,countEnh,'pr');plot(GE1,ydata22,'r:');
%title('GE vs DAR count');
 %title(['GE main vs DAR count, slope=',num2str(slope)]);

sl=slope;


%figure;
 %[muGE1,ce]=mean_Expr_CountEnh(countEnh,GE1);
 %[inter1,slope1,ydata1]=slope_int2(countEnh,GE1);
 % plot(countEnh,ydata1,'r');
 
 %[muGE2,cem]=mean_Expr_CountEnh2(countEnh,GE2);
 %[inter2,slope2,ydata2]=slope_int2(countEnh,GE2);
 %plot(countEnh,ydata2,'m');

 %title('countEnh: around genes with GE1 vs GE2');grid;
  
 
  
  end

end% function

%=====================subf 14 
function [inter,slope,ydata2]=slope_int2(x,y);


%------------------initiate output
inter=0;slope=0;ydata2=[];

xdata =x;% [1 2 3 4 5 6 7 8 9 10];
ydata =y;% [1.71 3.43 5.05 6.81 8.4 10.14 11.99 13.53 15.2 17.05];

%number of data points
mmax = length(x);
n=2;

  for i=1:n
    for j=1:n
    M(i,j) =0;
    end
    b(i) = 0;
  end

  for i=1:mmax
  M(1,1) = M(1,1) + 1;
  M(1,2) = M(1,2) + xdata(i);
  M(2,1) = M(2,1) + xdata(i);
  M(2,2) = M(2,2) + xdata(i)^2;
  b(1) = b(1) + ydata(i);
  b(2) = b(2) + ydata(i)*xdata(i);
  end

tolerance = 0.001;

max_it = 100;
a = [0 0];

 output = [0 a];

a=a';
b=b';

  for k=1:max_it

    anew(1) =  (b(1) - M(1,2:n)*a(2:n))/M(1,1);
    anew(n) = (b(n) - M(n,1:n-1)*anew(1:n-1)')/M(n,n);
    norm1 = max(abs(anew-a'));
    norm2 = max(abs(anew));
    if ((norm1/norm2)<tolerance)
        output = [output; k anew];
        break
    end
    output = [output; k anew];
    a=anew';
  end

iteration = k;

error=0;

  for i=1:mmax 
  ydata2(i) = a(1) + a(2)*xdata(i);
  error = error + (ydata(i) - ydata2(i))^2;
  end

%error
%figure;
%plot(xdata,ydata,'+');
%hold on;
%plot(xdata,ydata2,'r');

%result=a(2);
slope=a(2);
inter=a(1);
end %fun

%==================subfun 15
function [dist,me_sd_cv]=distances_analysis_np(Comm)
%3========================distances between genes

%------------------------outputs
%dist_p-distances without p-outpier (percentile)
%Yp distance giving p-percentile

%-----------------------inputs: thr_sig- just for plotting title

%Comm oor DiffEc etc
%%Chromosome	Start	End	Probe.Strand	E7.5_Ectoderm	E7.5_Endoderm	E7.5_Mesoderm
%1	3205901	3671498	0	04.35045	02.3982441	03.20762
%1	4343507	4360314	0	09.973417	08.495143	03.2783616

%--------------------intiate
me_sd_cv=[0,0,0];
dist=[];


 if min(size(Comm))>0,
starts1=Comm(:,2);
%chr=Comm(:,1);
dist=[];

%
 starts=sortrows(starts1);

  for i=1:length(starts)-1,
    dist(i)=abs(starts(i+1)-starts(i));
  end


%figure;histogram(dist_p,100);grid;
%title(['distribution of distances bw genes, thr=',num2str(thr_sig)]);

%-----------------mean sd without outlier
me=mean(dist);
sd=std(dist);

   if me~=0,
   cv=sd/me;
   else
    cv=-1;% bad
   end
  me_sd_cv=[me,sd,cv];

  end


end %func
%=================subf16
function [ng_LMH,ne_LMH,prop_LMH,TabH1,TabH2,TabH3,cGE1,cGE2,cGE3]=num_enh_lineage_11(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,endN,mesN,indEctN,indEndN,indMesN,vec_EcN,vec_EnN,vec_MN,dar1,dar2,dar3,chrN)
   
 % became -1  becaese of simon_11, because of ind_1 january 22, linea
% already separated by chrN though...

          %-------------------13-00 new
          [TabH1,TabH2,TabM1,TabM2,TabL1,TabL2,BDH1,BDH2,BDM1,BDM2,BDL1,BDL2,cGE1,cGE2,HE_bin_hist,ME_bin_hist,LE_bin_hist,re_rm,distH1,distH2,distM1,distM2,distL1,distL2,tar_darH,tar_darM,tar_darL,me12_HML,ncG1,ncG2,ngw1,ngw2,numLMH_12]=simon_plot_distHML_EcEn(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,endN,indEctN,indEndN,vec_EcN,vec_EnN,dar1,dar2,chrN);%,nameH,nameHN,nameL,nameLN);
          numEcEn_LMH=numLMH_12;

          [TabH1,TabH3,TabM1,TabM3,TabL1,TabL3,BDH1,BDH3,BDM1,BDM3,BDL1,BDL3,cGE1,cGE3,HE_bin_hist,ME_bin_hist,LE_bin_hist,re_rm,distH1,distH3,distM1,distM3,distL1,distL3,tar_darH,tar_darM,tar_darL,me13_HML,ncG1,ncG3,ngw1,ngw3,numLMH_13]=simon_plot_distHML_EcM(thr_vic,thr10,thrh,t_GE,te,gene_names,ectN,mesN,indEctN,indMesN,vec_EcN,vec_MN,dar1,dar3,chrN);%,nameH,nameHN,nameL,nameLN);
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
%==============subf17
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