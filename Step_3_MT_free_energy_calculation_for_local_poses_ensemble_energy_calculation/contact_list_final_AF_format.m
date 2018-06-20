
function torsion_list_num = contact_list_final_AF_format(torsion_list_final_AF)

%load torsion_list_final_AF
for i = 2:size(torsion_list_final_AF,1)
    torsion_name1 = torsion_list_final_AF{i,1};
    clear torsion_name_num1
    for j = 1:length(torsion_name1)
        if strcmp(torsion_name1(j),'C')
            torsion_name_num1(j) = '1';
        elseif strcmp(torsion_name1(j),'N')
            torsion_name_num1(j) = '2';
        elseif strcmp(torsion_name1(j),'O')
            torsion_name_num1(j) = '3';
        elseif strcmp(torsion_name1(j),'S')
            torsion_name_num1(j) = '4';
        elseif strcmp(torsion_name1(j),'A')
            torsion_name_num1(j) = '1';
        elseif strcmp(torsion_name1(j),'B')
            torsion_name_num1(j) = '2';
        elseif strcmp(torsion_name1(j),'D')
            torsion_name_num1(j) = '3';    
        elseif strcmp(torsion_name1(j),'E')
            torsion_name_num1(j) = '4';    
        elseif strcmp(torsion_name1(j),'G')
            torsion_name_num1(j) = '5';        
        elseif strcmp(torsion_name1(j),'H')
            torsion_name_num1(j) = '6';            
        elseif strcmp(torsion_name1(j),'Z')
            torsion_name_num1(j) = '7';    
        elseif strcmp(torsion_name1(j),'1')
            torsion_name_num1(j) = '1';    
        elseif strcmp(torsion_name1(j),'2')
            torsion_name_num1(j) = '2';
        elseif strcmp(torsion_name1(j),'3')
            torsion_name_num1(j) = '3';
        end
    end
    torsion_name1 = torsion_list_final_AF{i,2};
    
    if strcmp(torsion_name1, 'ALA')==1
        torsion_res_num1 = 101;
    elseif strcmp(torsion_name1, 'ARG')==1
        torsion_res_num1 = 102;
    elseif strcmp(torsion_name1, 'ASN')==1 || strcmp(torsion_name1, 'ASX')==1
        torsion_res_num1 = 103;
    elseif strcmp(torsion_name1, 'ASP')==1
        torsion_res_num1 = 104;
    elseif strcmp(torsion_name1, 'CYS')==1
        torsion_res_num1 = 105;
    elseif strcmp(torsion_name1, 'GLN')==1
        torsion_res_num1 = 106;
    elseif strcmp(torsion_name1, 'GLU')==1 || strcmp(torsion_name1, 'GLX')==1
        torsion_res_num1 = 107;
    elseif strcmp(torsion_name1, 'GLY')==1
        torsion_res_num1 = 108;
    elseif strcmp(torsion_name1, 'HIS')==1 || strcmp(torsion_name1, 'HID')==1 || strcmp(torsion_name1, 'HIE')==1 || strcmp(torsion_name1, 'HIP')==1
        torsion_res_num1 = 109;
    elseif strcmp(torsion_name1, 'ILE')==1
        torsion_res_num1 = 110;
    elseif strcmp(torsion_name1, 'LEU')==1
        torsion_res_num1 = 111;
    elseif strcmp(torsion_name1, 'LYS')==1
        torsion_res_num1 = 112;
    elseif strcmp(torsion_name1, 'MET')==1
        torsion_res_num1 = 113;
    elseif strcmp(torsion_name1, 'PHE')==1
        torsion_res_num1 = 114;
    elseif strcmp(torsion_name1, 'PRO')==1
        torsion_res_num1 = 115;
    elseif strcmp(torsion_name1, 'SER')==1
        torsion_res_num1 = 116;
    elseif strcmp(torsion_name1, 'THR')==1
        torsion_res_num1 = 117;
    elseif strcmp(torsion_name1, 'TRP')==1
        torsion_res_num1 = 118;
    elseif strcmp(torsion_name1, 'TYR')==1
        torsion_res_num1 = 119;
    elseif strcmp(torsion_name1, 'VAL')==1
        torsion_res_num1 = 120;
    else
        torsion_res_num1 = 121;
    end
    

    
    
    
    
    
    
    
    
    
    
    torsion_name2 = torsion_list_final_AF{i,3};
    clear torsion_name_num2
    for j = 1:length(torsion_name2)
        if strcmp(torsion_name2(j),'C')
            torsion_name_num2(j) = '1';
        elseif strcmp(torsion_name2(j),'N')
            torsion_name_num2(j) = '2';
        elseif strcmp(torsion_name2(j),'O')
            torsion_name_num2(j) = '3';
        elseif strcmp(torsion_name2(j),'S')
            torsion_name_num2(j) = '4';
        elseif strcmp(torsion_name2(j),'A')
            torsion_name_num2(j) = '1';
        elseif strcmp(torsion_name2(j),'B')
            torsion_name_num2(j) = '2';
        elseif strcmp(torsion_name2(j),'D')
            torsion_name_num2(j) = '3';    
        elseif strcmp(torsion_name2(j),'E')
            torsion_name_num2(j) = '4';    
        elseif strcmp(torsion_name2(j),'G')
            torsion_name_num2(j) = '5';        
        elseif strcmp(torsion_name2(j),'H')
            torsion_name_num2(j) = '6';            
        elseif strcmp(torsion_name2(j),'Z')
            torsion_name_num2(j) = '7';    
        elseif strcmp(torsion_name2(j),'1')
            torsion_name_num2(j) = '1';    
        elseif strcmp(torsion_name2(j),'2')
            torsion_name_num2(j) = '2';
        elseif strcmp(torsion_name2(j),'3')
            torsion_name_num2(j) = '3';
        end
    end
    
    
    torsion_name2 = torsion_list_final_AF{i,4};
    
    if strcmp(torsion_name2, 'ALA')==1
        torsion_res_num2 = 101;
    elseif strcmp(torsion_name2, 'ARG')==1
        torsion_res_num2 = 102;
    elseif strcmp(torsion_name2, 'ASN')==1 || strcmp(torsion_name2, 'ASX')==1
        torsion_res_num2 = 103;
    elseif strcmp(torsion_name2, 'ASP')==1
        torsion_res_num2 = 104;
    elseif strcmp(torsion_name2, 'CYS')==1
        torsion_res_num2 = 105;
    elseif strcmp(torsion_name2, 'GLN')==1
        torsion_res_num2 = 106;
    elseif strcmp(torsion_name2, 'GLU')==1 || strcmp(torsion_name2, 'GLX')==1
        torsion_res_num2 = 107;
    elseif strcmp(torsion_name2, 'GLY')==1
        torsion_res_num2 = 108;
    elseif strcmp(torsion_name2, 'HIS')==1 || strcmp(torsion_name2, 'HID')==1 || strcmp(torsion_name2, 'HIE')==1 || strcmp(torsion_name2, 'HIP')==1
        torsion_res_num2 = 109;
    elseif strcmp(torsion_name2, 'ILE')==1
        torsion_res_num2 = 110;
    elseif strcmp(torsion_name2, 'LEU')==1
        torsion_res_num2 = 111;
    elseif strcmp(torsion_name2, 'LYS')==1
        torsion_res_num2 = 112;
    elseif strcmp(torsion_name2, 'MET')==1
        torsion_res_num2 = 113;
    elseif strcmp(torsion_name2, 'PHE')==1
        torsion_res_num2 = 114;
    elseif strcmp(torsion_name2, 'PRO')==1
        torsion_res_num2 = 115;
    elseif strcmp(torsion_name2, 'SER')==1
        torsion_res_num2 = 116;
    elseif strcmp(torsion_name2, 'THR')==1
        torsion_res_num2 = 117;
    elseif strcmp(torsion_name2, 'TRP')==1
        torsion_res_num2 = 118;
    elseif strcmp(torsion_name2, 'TYR')==1
        torsion_res_num2 = 119;
    elseif strcmp(torsion_name2, 'VAL')==1
        torsion_res_num2 = 120;
    else
        torsion_res_num2 = 121;
    end
    
    torsion_list_num(i,1) = str2num(torsion_name_num1);
    torsion_list_num(i,2) = torsion_res_num1;
    torsion_list_num(i,3) = str2num(torsion_name_num2);
    torsion_list_num(i,4) = torsion_res_num2;
end
            
            
            
            
            
            
    