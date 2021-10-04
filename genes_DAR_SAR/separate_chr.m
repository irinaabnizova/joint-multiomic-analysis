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
