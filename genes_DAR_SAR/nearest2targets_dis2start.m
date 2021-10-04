function [count_Enh_target,EnhTarget,TargetEnh,dis2start,perc_body,ind_Enh]=nearest2targets_dis2start(thr_vic,Enh,target,add)

%=======compute  occurence of DARs within visinity of target TSS
%---add=shift in bp from start-end of target, if needed
%targets=genes/enhancers/TSS etc

%-------------------Input bed files: or more columns: chr_num start_pos end_pos

%7  3097633 200 0.44 0.15 36 93 6.000000e-01 2611
%7   3106133 200 0.25 0.10 72 84 4.600000e-01 2044
%7   3149333 400 0.35 0.23 40 31 1.740000e+00 2375

%--------------------output
%count_Enh_target= count of enhancers around each target pos: vector of
%length (=target max size)

%---------------------initiate
count_Enh_target=[];
EnhTarget=[];
TargetEnh=[];

dis2start=[];
perc_body=[];
ind_Enh=[];

%ind_targetEnh=[];
  
 
%-----------------------------if non-empty target and Enh data
  
if max(size(target)) & max(size(Enh)) >0,

    % -----------------------------------pos=starts
    pos_target=target(:,2);
    posEnh=Enh(:,2);
    end_target=target(:,3);

    
    %--initiate number of enhancers around each target start position
    for j=1:length(pos_target),
        count_Enh_target(j)=0;
    end
    
 
    cenh=0;% count event of target having Enh around
    
    % count all enhancers around targets: all sorted!
 
    for j=1:length(pos_target),
       % for sorted lists 
       for i=j:length(posEnh),
           % if Enhancer is close enough to target(j) start
            if abs(posEnh(i)-pos_target(j)) < thr_vic,
                count_Enh_target(j)=count_Enh_target(j)+1;% per target: count enhancers located in its thr-neighbourhood
                cenh=cenh+1;% count events when enhancers located in the thr-neibourhood of chrN target targets
                dis2start(cenh)=posEnh(i)-pos_target(j);
                len_target(cenh)=-pos_target(j)+end_target(j);% for inside normalization
  
                %Enh=DARs here within targets vicinity
                EnhTarget(cenh,:)=Enh(i,:);% set of Enhancer around target
                TargetEnh(cenh,:)=target(j,:);
                ind_Enh(cenh)=i;% index in the Enh file
            end
       end % i
   end % j


  %----------------------percentage pos of region within target body +-add
  nc=0;
  perc_body=[];
  %add=600;
  for i=1:length(len_target),
      if dis2start(i)>= (0 -add ) &&  dis2start(i)<= (len_target(i)+add),
          nc=nc+1;
          perc_body(nc)= round(100*dis2start(i)/len_target(i));% pos/length
      end
  end


end % main if non empty data target and Enh
