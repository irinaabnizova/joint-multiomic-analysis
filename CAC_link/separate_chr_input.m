function [dar1,dar2,dar3,ectN,indEctN,vec_EcN,endN,indEndN,vec_EnN,mesN,indMesN,vec_MN,ng123,ne123] = separate_chr_input( ectMZ,endMZ,mesMZ,DAR1,DAR2,DAR3,chrN )


%--------------------outputs
    %ng123=numGenes_EcEnM_chrN;
    %ne123=num_DAR_EcEnM_chrN;
    ng123=0;ne123=0;
    
    [dar1,dar3]=dars_chr_2(DAR1,DAR3, chrN);
    [dar1,dar2]=dars_chr_2(DAR1,DAR2, chrN);
    
    [ectN,indEctN]=separate_chr_ind(ectMZ,ectMZ(:,8), chrN);
    [endN,indEndN]=separate_chr_ind(endMZ,endMZ(:,8), chrN);
    [mesN,indMesN]=separate_chr_ind(mesMZ,mesMZ(:,8), chrN);

    
    vec_EcN=ectN(:,5);
    vec_EnN=endN(:,6);
    vec_MN=mesN(:,7);

    numGenes_EcEnM_chrN=[length(vec_EcN),length(vec_EnN),length(vec_MN)];

    %-------------------check for non-emptiness !!!
    num_DAR_EcEnM_chrN=[length(dar1(:,1)),length(dar2(:,1)),length(dar3(:,1))];

    ng123=numGenes_EcEnM_chrN;
    ne123=num_DAR_EcEnM_chrN;
end
%===================subfunc
function [dare,darc,indE,indC]=dars_chr_2(DARE4_100,DARC4_100, chrN)

%2.1-------------------two dars for one chrom
dare=[];darc=[];indE=[];indC=[];
[DARE4_100_chr,indE]=separate_chr(DARE4_100, chrN);
[DARC4_100_chr,indC]=separate_chr(DARC4_100, chrN);

  if size(DARE4_100_chr)>0 & size(DARE4_100_chr)>0,
  k=min(length(DARE4_100_chr(1,:)),length(DARC4_100_chr(1,:)));% number of columns
    %-------------only positions: start end score
  dare=DARE4_100_chr(:,2:k);
  darc=DARC4_100_chr(:,2:k);
  end
end

%======subfunc
function [chr7_coor,ind7_enh]=separate_chr(chr_enh, chrNum)
 %separate one numbered chromosome for enhancer table data, keeping indexes
 %in that table

%-------------------Input
%chr_gene
% chr gene_start  )
%1	  3062534	
%1	  3071886	
%....
%666(Y)	3191164	
%66 (X)	3399979	

%-----------------Output
%---------------two/three column matrix, if chrN=chr7
%chr7_coor= 
% chr gene_start  )
%7	  3062534	
%7	  3071886	
%....
%7	3191164	
%7	3399979	

%------------------------------body

 chr7_coor=[];%extracted chr7 : 2 columns
%----------------vector
 ind7_enh=[];% their indexes in the original table 

 c7=0;

  for i=1:length(chr_enh(:,1)),
    if chr_enh(i,1)==chrNum,
        c7=c7+1;
        chr7_coor(c7,:)=chr_enh(i,:);
        ind7_enh(c7)=i;%keep this index from the original all-chrom file
    end
  end
%c7% 20
end
%======================subfunc
function [chrN_coor,indN_enh]=separate_chr_ind(chr_enh,ind_enh, chrNum)

%separate one numbered chromosome for enhancer table data and enhancer index!!!!, keeping indexes
%in that table

%-------------------Input
%chr_gene
% chr gene_start  )
%1	  3062534	
%1	  3071886	
%....
%666(Y)	3191164	
%66 (X)	3399979	

%-----------------Output
%---------------two/three column matrix, if chrN=chr7
%chr7_coor= 
% chr gene_start  )
%7	  3062534	
%7	  3071886	
%....
%7	3191164	
%7	3399979	

%------------------------------body

chrN_coor=[];%extracted chrN : 2 columns
%----------------vector
indN_enh=[];% their indexes in the given table 

cN=0;

  for i=1:length(chr_enh(:,1)),
    if chr_enh(i,1)==chrNum,
        cN=cN+1;
        chrN_coor(cN,:)=chr_enh(i,:);
        indN_enh(cN)=ind_enh(i);%keep this index from the already Dif separated, all-chrom data
    end
  end
%cN% 20
end