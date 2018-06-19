function [ligand_final_str,ligand_final_dG,ligand_final_repulse_num] = ErelaxB(ligand_select_X1,dG_select_pass_X1,index,pocket,vdw_d)

ligand_select_coor_X1 = ligand_select_X1(dG_select_pass_X1(index,3)).coordinate;
ligand_select_cent_X1(1,1) = mean(ligand_select_coor_X1(:,1));
ligand_select_cent_X1(1,2) = mean(ligand_select_coor_X1(:,2));
ligand_select_cent_X1(1,3) = mean(ligand_select_coor_X1(:,3));

lig_sel_modi(:,1) = ligand_select_coor_X1(:,1)+ 0-ligand_select_cent_X1(1,1);
lig_sel_modi(:,2) = ligand_select_coor_X1(:,2)+ 0-ligand_select_cent_X1(1,2);
lig_sel_modi(:,3) = ligand_select_coor_X1(:,3)+ 0-ligand_select_cent_X1(1,3);
lig_sel_modi(:,4:10) = ligand_select_coor_X1(:,4:10);

pocket_new(:,1) = pocket(:,1)+ 0-ligand_select_cent_X1(1,1);
pocket_new(:,2) = pocket(:,2)+ 0-ligand_select_cent_X1(1,2);
pocket_new(:,3) = pocket(:,3)+ 0-ligand_select_cent_X1(1,3);
pocket_new(:,4:10) = pocket(:,4:10);

ligand_select_repulse_num_X1 = ligand_select_X1(dG_select_pass_X1(index,3)).repulse;

[ligand_new_str_X1] = strEsearch_multiV2(lig_sel_modi,pocket_new,dG_select_pass_X1(index,1),ligand_select_repulse_num_X1,vdw_d);
clear dG_select_new
for i2 = 1:size(ligand_new_str_X1,2)
    if ~isempty(ligand_new_str_X1(i2).dG)
        dG_select_new_X1(i2,1)=ligand_new_str_X1(i2).dG(1,1);
    elseif isempty(ligand_new_str_X1(i2).dG)
        dG_select_new_X1(i2,1)=0;
    end
end
dG_select_new_X1(:,2)=1:size(dG_select_new_X1,1);
dG_select_new_X1 = sortrows(dG_select_new_X1,[1]);
dG_select_pass_new_X1=dG_select_new_X1(1,:);
ligand_select_repulse_num_new_X1 = ligand_new_str_X1(dG_select_pass_new_X1(1,2)).repulse_num;
lig_sel_modi_new_X1 = ligand_new_str_X1(dG_select_pass_new_X1(1,2)).structure;
[ligand_final_str,ligand_final_dG,ligand_final_repulse_num,loop] = strEsearch(lig_sel_modi_new_X1,pocket_new,dG_select_pass_new_X1(1,1),ligand_select_repulse_num_new_X1,vdw_d);

ligand_final_str(:,1) = ligand_final_str(:,1) + ligand_select_cent_X1(1,1);
ligand_final_str(:,2) = ligand_final_str(:,2) + ligand_select_cent_X1(1,2);
ligand_final_str(:,3) = ligand_final_str(:,3) + ligand_select_cent_X1(1,3);
ligand_final_str(:,4:10)=ligand_final_str(:,4:10);


