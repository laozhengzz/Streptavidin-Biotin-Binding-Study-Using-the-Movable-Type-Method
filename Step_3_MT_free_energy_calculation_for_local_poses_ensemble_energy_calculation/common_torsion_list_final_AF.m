function  [common_list1,common_list1_Namelist,common_list2,common_list2_Namelist] = common_torsion_list_final_AF(list1,list1_Namelist,list2,list2_Namelist)

[C,ia,ib] = intersect(list1(:,10:11),list2(:,10:11),'rows');

common_list1 = list1(ia,:);
common_list2 = list2(ib,:);
common_list1_Namelist = list1_Namelist(ia,:);
common_list2_Namelist = list2_Namelist(ib,:);

















