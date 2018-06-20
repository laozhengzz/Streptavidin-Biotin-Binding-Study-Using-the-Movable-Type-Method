%function [P_lowest, PL_lowest P_Ensemble, PL_Ensemble] = streptavidin_1dpmf_boltzenergy(protein_dir) 
clear all;
fp=fopen(strcat('P_PL_1d_pmf.txt'),'wt');


for a=1:819
clear E_P_CL_ind E_PL_CL_ind;
protein_dir=sprintf('/mnt/home/bansalnu/Protein_ligand_energy_code/cluster_%d/data.txt',a)
load(protein_dir);
RcutoffPP = 6;
RcutoffPL = 6;
%P_Tor	P_Non	ZE_P_san	ZE_L_san	ZE_P_water	ZE_L_water	COM 	ZE_PL_san	ZE_PL_water

E_P_CL_ind(:,1) = data(:,1)+data(:,2)./(RcutoffPP*10) + (data(:,3) + data(:,4));
E_PL_CL_ind(:,1) = data(:,1)+data(:,2)./(RcutoffPP*10) + data(:,7)./(RcutoffPL)  + data(:,8) + data(:,9)*-0.05;


temp_list1 = sortrows(E_P_CL_ind(:,1),[1]);
E_P_CL_ind_low(:,1) = temp_list1(1:10,:)-mean(temp_list1(1:10,:));

temp_list2 = sortrows(E_PL_CL_ind(:,1),[1]);
E_PL_CL_ind_low(:,1) = temp_list2(1:10,:)-mean(temp_list2(1:10,:));


Ensemble_E(1,1)=sum(temp_list1(1:10,:).*exp(E_P_CL_ind_low(:,1)./-0.5918)./sum(exp(E_P_CL_ind_low(:,1)./-0.5918)));%+mean(temp_list1(1:100,:));
Ensemble_E(2,1)=sum(temp_list2(1:10,:).*exp(E_PL_CL_ind_low(:,1)./-0.5918)./sum(exp(E_PL_CL_ind_low(:,1)./-0.5918)));%+mean(temp_list2(1:100,:));

State_E_min(1,1) = min( E_P_CL_ind(:,1));
State_E_min(2,1) = min( E_PL_CL_ind(:,1));


[sem1,r1] = find(E_P_CL_ind(:,1) == min( E_P_CL_ind(:,1)));
[sem2,r2] = find(E_PL_CL_ind(:,1) == min( E_PL_CL_ind(:,1)));

P_lowest(a,1) = State_E_min(1,1);
PL_lowest(a,1) = State_E_min(2,1);

P_Ensemble(a,1) = Ensemble_E(1,1);
PL_Ensemble(a,1) = Ensemble_E(2,1);
fprintf(fp,' % 5d  % 5d  % 5d  % 5d\n',P_lowest(a,1), PL_lowest(a,1), P_Ensemble(a,1), PL_Ensemble(a,1));

end

