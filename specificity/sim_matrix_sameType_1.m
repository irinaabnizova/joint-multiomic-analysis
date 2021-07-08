function [Mpp]=sim_matrix_sameType_1(fp_listSEGi,fp_listEcti,fp_listEndi,fp_listMesi)
%----------------------Compute sililarity matrix for promoters
%display(['factor=',name2str(fact)])

feature4ip=fp_listSEGi;

feature1ip=fp_listEcti;
feature2ip=fp_listEndi;
feature3ip=fp_listMesi;


%display('similarity between SEGs and DEGs sequences: Ect End Mes');
RMSE12pp = sqrt(mean((feature1ip-feature2ip).^2));
RMSE13pp = sqrt(mean((feature1ip-feature3ip).^2));
RMSE14pp = sqrt(mean((feature1ip-feature4ip).^2));

%display('diagonal: similarity between SEGs and DEGs sequences themselves');

RMSE11pp = sqrt(mean((feature1ip-feature1ip).^2));
RMSE22pp = sqrt(mean((feature2ip-feature2ip).^2));
RMSE33pp = sqrt(mean((feature3ip-feature3ip).^2));
RMSE44pp = sqrt(mean((feature4ip-feature4ip).^2));

%display('similarity between DEGs sequences themselves: Ect End Mes');

RMSE23pp = sqrt(mean((feature2ip-feature3ip).^2));
RMSE24pp = sqrt(mean((feature2ip-feature4ip).^2));
RMSE34pp = sqrt(mean((feature3ip-feature4ip).^2));

Mpe(1,1)=RMSE11pp;
Mpp(1,2)=RMSE12pp;
Mpp(1,3)=RMSE13pp;
Mpp(1,4)=RMSE14pp;

Mpp(2,1)=RMSE12pp;
Mpp(3,1)=RMSE13pp;
Mpp(4,1)=RMSE14pp;


Mpp(2,2)=RMSE22pp;
Mpp(2,3)=RMSE23pp;
Mpp(2,4)=RMSE24pp;

Mpp(3,2)=RMSE23pp;
Mpp(4,2)=RMSE24pp;


Mpp(3,3)=RMSE33pp;
Mpp(3,4)=RMSE34pp;
Mpp(4,4)=RMSE44pp;

Mpp(4,3)=RMSE34pp;

%Mpp
