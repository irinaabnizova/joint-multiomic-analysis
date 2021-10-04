function [bin, count_bin,Ncb,Ncbv]=histo5K_normVicNov(thr_vic,d2s1,num_genes,ngw)

% -9 Nov 2020- precompute hist (Ncb,Ncbv) as 0
 %compute  counts (count_bin) of DAR starts in 1KB up/downstream TSSs from d3s 
 %Ncb= Normalised count bin

 % ngw total number of genes, includingthose in the vicinity of TSS: it
 % it might be >num_genes
 
 %Ncbv=count_bin/ngw; Normalised count bin vicinity
 
 bin=-thr_vic:5000:thr_vic;
 
 %------initialiise per bin counts
 
 for j=1:length(bin),
     count_bin(j)=0;
     Ncb(j)=0;
     Ncbv(j)=0;
 end
 
 
 if size(d2s1)>0 & num_genes>0 & ngw>0,
 
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
 end
