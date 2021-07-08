function [matrix_sim_Original,vectors_sim_Shuffled,specificity,prob]=MAIN_compute_specificity(FLCC,list_SEG,list_Ect,list_End,list_Mes,fS,quart_weights)
%---------------------
%Computes original and permuted (100 time) matrix of promoter-promoter
%similarities within DEGs and between DEGs SEGs

%-----------------------INPUTS- names of files!!!!
%FLC=Dictionary of all possible here motifs

%-----------------------DEfault option
if nargin < 7,
  display('Default Parameters of specificity computing from four lists of TF motifs and motif Dictionary(FLC)');
  fS=0.5
  %--------------------linear
  t=0.2;
  va1=8*t;
  va2=6*t;
  va3=2*t;
  va4=0.1*t;
  quart_weights=[va1,va2,va3,va4]
end

%-----------------------open file names

    %[FLCC,list_SEG,list_Ect,list_End,list_Mes,outOriginal,outShuffled]=test_enhancer_names()
   
    %=========================check if input data is OK, and load it
    ew=fopen(FLCC);
    if ew<0,
        display('wrong/non existent input1 Dictionary (FLC) file');
        return
    end
    fclose(ew);
    FLC=importdata(FLCC);%alphabetical, 277
   
    %---------------------------------SEG regulatory regions
    e_seg=fopen(list_SEG);
    if e_seg<0,
        display('wrong/non existent input2 seg list file');
        return
    end
    fclose(e_seg);
    listS1=importdata(list_SEG);


    %---------------------------------DEG regulatory regions
    e_ect=fopen(list_Ect);
    if e_ect<0,
        display('wrong/non existent input3 ect list file');
        return
    end
    fclose(e_ect);
    list11=importdata(list_Ect);

    e_end=fopen(list_End);
    if e_end<0,
        display('wrong/non existent input4 end list file');
        return
    end
    fclose(e_end);
    list21=importdata(list_End);
    

    e_mes=fopen(list_Mes);
    if e_mes<0,
        display('wrong/non existent input5 mes list file');
        return
    end
    fclose(e_mes);
    list31=importdata(list_Mes);
    


%--------------------main body
ind=1:length(FLC);


%---------------------control list length
f1=fS;
f2=fS;
f3=fS;
listS=listS1(1:round(fS*length(listS1)));
indS=1:length(listS);
list1=list11(1:round(f1*length(list11)));
ind1=1:length(list1);
list2=list21(1:round(f2*length(list21)));
ind2=1:length(list2);
list3=list31(1:round(f3*length(list31)));
ind3=1:length(list3);

numSEcEnM_controlLength=[length(indS),length(ind1),length(ind2),length(ind3)]

%---------------use quartile weights: usually linear or exp
%with first quatrile giving max weight

%quart_weights=[va1,va2,va3,va4]
va1=quart_weights(1);
va2=quart_weights(2);
va3=quart_weights(3);
va4=quart_weights(4);

%---------------------compute permuted and original feature-vectors
%for each list
N=100
[M2pp,rel_err,DSP,errP,pr]=permute_original_loop(N,va1, va2, va3,va4,FLC,ind,listS,indS,list1,ind1,list2,ind2,list3,ind3);
%DSP
rel_err_plot=ones(N,1)*rel_err;
figure;histogram(errP,30);title('distribution of permuted DD-DS scores differences, DARs enhancers');
hold on;
histogram([errP rel_err_plot'],30);
histogram(errP,30);
xlabel('score-differences');

%---------------------compute number of errP >= rel_err
nm=0;
for i=1:N,
    if errP(i)>=rel_err,
        nm=nm+1;
    end
end
prob=nm/N
me_permutedDSP=mean(DSP)
matrix_sim_Original=M2pp;
vectors_sim_Shuffled=DSP;
specificity=rel_err

%===============write Original into output files, with accessibility level per lineage
   [a]=save_chromatin_layers(DAR,outDAR);
   [b]=save_chromatin_layers(SARH,outSAR);
 