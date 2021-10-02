function [count_Enh_target,EnhTarget,TargetEnh,num_int_Enh_tar,per_int_Enh_tar,dis2start,perc_body,ind_Enh]=around_target_dis2start(thr_vic,Enh,target,add);

% to compute relatitive inside target location of enhancers: 

%dis2start= only for intersected enhs around target

%---------------------initiate
count_Enh_target=[];
EnhTarget=[];
TargetEnh=[];
num_int_Enh_tar=[];
per_int_Enh_tar=[];
dis2start=[];
perc_body=[];
ind_Enh=[];

%--------------
if size(Enh)>0 & size(target) >0, 
    num_tar=length(target(:,1));
    num_Enh=length(Enh(:,1));


   [count_Enh_target,EnhTarget,TargetEnh,dis2start,perc_body,ind_Enh]=nearest2targets_dis2start(thr_vic,Enh,target,add);


   %EnhTarget,TargetEnh- interaction matrices-types

 
   %--------------------number of not covered target genes
 
    nz=0;z=0;
    count_Enh_targetNZ=[];
    for i=1:length(count_Enh_target),
        if count_Enh_target(i)==0,
           nz=nz+1;
        else
           z=z+1;
           count_Enh_targetNZ(z)=count_Enh_target(i);
        end
    end
    num_tar_int=z;% non empty interacting targets


   per_tar=100*(length(count_Enh_target)-nz)/length(target(:,1));


  %  int= interaction here:  num_Enh_int = number enhancers interacting
  %  with target: likely duplicated sometimes (the more vicinity!)
   if max(size(EnhTarget)) > 0,
      num_Enh_int=length(EnhTarget(:,1));
   else
      num_Enh_int=0;
   end
 
   per_Enh=100*num_Enh_int/num_Enh;
   input_Enh_tar=[num_Enh,num_tar];
   num_int_Enh_tar=[num_Enh_int,num_tar_int];
   per_int_Enh_tar=round([per_Enh,per_tar]);
  
end % non empty data