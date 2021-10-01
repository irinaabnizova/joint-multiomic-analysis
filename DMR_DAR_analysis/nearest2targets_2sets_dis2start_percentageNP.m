function [count_Enh_target,EnhTarget,TargetEnh,dis2start,perc,ind_Enh]=nearest2targets_2sets_dis2start_percentageNP(thr_vic,Enh,target,add)

%=======compute relative occurence of DARs within K4me3 bodies, =from the strats!!!!
% not swapped start and ends in targets

%---add=shift in bp from start-end of target, if needed

%=======find nearest ROI in  thr=1000 bp (?) of target, only positions (no
%indexes)

%targets=genes/enhancers/TSS etc
%etc
%-------------------Input target=3 or more columns: chr_num start_pos end_pos

 %  with added chr first column!!!!!
 %7  3097633 200 0.44 0.15 36 93 6.000000e-01 2611
%7   3106133 200 0.25 0.10 72 84 4.600000e-01 2044
%7   3149333 400 0.35 0.23 40 31 1.740000e+00 2375
%3160133 400 0.16 0.45 31 29 -2.940000e+00 -3442
%3173833 300 0.25 0.00 28 22 9.300000e-01 2583


%Enh=starts of DARs/enhancers  Signif windows starts: samr format! chr(!!!) positions


%--------------------output
%count_Enh_target= count of enhancers around each target pos: vector of target length
%pos_EnhTarget



%==========================from the starts
%K4me3=load('C:\Users\ia1\Documents\MATLAB\molec_layers\analysis_joint\intersect_correlate\data_gene_enh\postions_K4me3_common_num.txt');
 %chrN=2;
%[K4me3_chr,ind_K4me3]=separate_chr(K4me3,chrN);%pos_enhM,pos_enhEc,pos_enhEn);

%targetE=K4me3_chr;
%targetM=K4me3_chr;

%---------------------initiate
count_Enh_target=[];
EnhTarget=[];
TargetEnh=[];
num_int_Enh_tar=[];
per_int_Enh_tar=[];
dis2start=[];
perc=[];
ind_Enh=[];

TargetEnh=[];
EnhTarget=[];
perc=[];
 pos_target=[];
  end_target=[];
   posEnh=[];
   len=[];

   ind_targetEnh=[];
% if target was start
  pos_target=target(:,2);
  posEnh=Enh(:,2);
  end_target=target(:,3);
  
  count_Enh_target=[];
  dis2start=[];% distances for enhancers relative to nearest gene start
  pos_tar_enh=[];
  
  
%-----------------------------if non-empty data
  
if max(size(target))>0,
  
    
    %--initiate number of enhancers around each target start position
    for j=1:length(pos_target),
      count_Enh_target(j)=0;
    end
    
 
    cenh=0;% count event of target having Enh around
    
    % count all enhancers around targets
    %ind_Enh=[];
    %ind_targetEnh=[];
   
   for j=1:length(pos_target),
    
       for i=1:length(posEnh),
          if abs(posEnh(i)-pos_target(j)) < thr_vic,
            count_Enh_target(j)=count_Enh_target(j)+1;% per target: count enhancers located in its thr-neighbourhood
            cenh=cenh+1;% count events when enhancers located in the thr-neibourhood of chrN target targets
            dis2start(cenh)=posEnh(i)-pos_target(j);
            len(cenh)=-pos_target(j)+end_target(j);% for inside normalization
  
            % =new insert 29 July: inside target position
              %   if dis2start(cenh)>= pos_target(j) && dis2start(cenh)<= end_target(j),%covered
              %   nc=nc+1;
               %  perc(nc)=100*dis2start(cenh)/len(cenh);% pos/length
               %  end
             
   
   %======================
            %Enh=DARs here within targets vicinity
            EnhTarget(cenh,:)=Enh(i,:);% set of Enhancer around target
            TargetEnh(cenh,:)=target(j,:);
            ind_Enh(cenh)=i;% index in the Enh
            %ind_targetEnh(cenh)=ind_target(j);% Enh-related indexes in the file of targets
         end
       end % i
   end % j

%---if swappwd strats ends: looks inside the body from the end point


%----------------------percentage pos of region within target body
nc=0;
perc=[];
%add=600;
for i=1:length(len),
    if dis2start(i)>= (0 -add ) &&  dis2start(i)<= (len(i)+add),
        nc=nc+1;
        perc(nc)= round(100*dis2start(i)/len(i));% pos/length
    end
end
%nc

    
%figure;hist(perc,20);title('distribution of occurences of DARs within body, relative');
%xlabel('percentage of  K4me3 body');
end % main if non empty data target
