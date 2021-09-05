function [pu,ch,amp]=get_pu_ch_meth_3(tot_access)

%compute Pu and Chas for each non empty access window

%-with ssc in case all elements are zero

%tot_access=   
%tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me'  BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%  tot=[winm' levBE' levBEn' levBC' covBE' covBEn' covBC' BE_npos_me' BE_npos_unme' BEn_npos_me'  BEn_npos_unme' BC_npos_me' BC_npos_unme' indm'];
%        1     2      3      4        5      6       7       8          9
%                                                                                      10          11              12            13     14 
    %5=====================compute Chastity/Purity of tot(:,2:4) (three
    %access Levels
  

%  tot_accesss uniquely diff expressed (as higher than other layers)
%  lineage-specific, same time e.g. E7.5

%--med = median, not mean
%tot_access_expr=load('C:\Users\ia1\Documents\MATLAB\molec_layers\analysis_joint\intersect_correlate\data_tot_access_enh\e75_expressions_chr_pos.txt');
%Chromosome	Start	End	Probe.Strand	E7.5_Ectoderm	E7.5_Endoderm	E7.5_Mesoderm
%1	3205901	3671498	0	04.35045	02.3982441	03.20762
%1	4343507	4360314	0	09.973417	08.495143	03.2783616


pu=[];
ch=[];
amp=[];

pup=[];
chp=[];

%meds=[median(pu_am_ec),median(pu_am_en),median(pu_am_m)]

%----------------desinged for three layers currently, e.g. E7.5 ect end mes

%============================INPUTS
%tot_access expression matrix with cols 5,6,7 are ect end mes
%ecto,endo,meso---------vectors of gE

%[tot_accessEc,tot_accessEn,tot_accessM,tot_access]=GE_lineage_EcEnM_sw(tot_access_expr);
%tot_access=Comm;

ssc=0.001;

ecto=1-tot_access(:,2);
endo=1-tot_access(:,3);
meso=1-tot_access(:,4);

%1----------------------the first round, to learn real amplitude
layers=[ecto endo meso];%tot_access(:,2:4);
nn=0;
for i=1:length(ecto),
    vec=layers(i,:);
    mi=min(vec);
    if mi > 0,
    pu(i)=max(vec)/sum(vec);
    amp(i)=max(vec);
    ch(i)=max(vec)/(max(vec)+second_max(vec));
    else
        % use ssc in case all vec elements are zero!
        nn=nn+1;
        con=abs(mi)+ssc;
        amp(i)=max(vec)+ssc;
        vec1=vec+con;
        if sum(vec1)>0,
        pu(i)=max(vec1)/sum(vec1);
        ch(i)=max(vec1)/(max(vec1)+second_max(vec1));
        else
        pu(i)=0;
        ch(i)=0;
        end
    end
end
    
%----------------------------round two-- DO NOT NEED for access???
%layers=tot_access(:,5:7);
%nn=0;
%for i=1:length(ecto),
 %   vec=layers(i,:);
  %  vecp=vec+(12-amp(i))^2;% to keep sharp: the more amp, the less is addition
    %vecp=layers_pos(i,:);
   % pup(i)=max(vecp)/sum(vecp);
   % chp(i)=max(vecp)/(max(vecp)+second_max(vecp));
%end

figure;
subplot(2,1,1);
histogram(pu);
title('pu across three layers');
xlabel('Pu of GE across Ect End Mes, each access win');
xlim([0 1]);

   subplot(2,1,2);
histogram(ch);
title('chastity across three layers');
 xlim([0 1]);
 
 


%==========================get sets of Diff tot_accesss
%thrPu=0.6;
%thramp=0.2;
