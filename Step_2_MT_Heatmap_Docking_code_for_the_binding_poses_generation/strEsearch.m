function [ligand_new_str_structure,ligand_new_str_dG,ligand_new_str_repulse_num,loop] = strEsearch(ligand_input,pocket_new,dG_select_pass,repulse_num_loop,vdw_d)
E_evaluate=0;
degree = -pi/18:pi/18:pi/18;
RcutoffPL = 6;
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;

const=halva*2;

loop=1;
while loop<20 
    
    
    E_evaluate=dG_select_pass;
    ln=1;
    clear ligand_new_all
    for axis1 = 1:3
        for dice = -0.1:0.1:0.1
            structure_tran = ligand_input;
            structure_tran(:,axis1)=ligand_input(:,axis1)+dice;
            for axis2 = 1:3
                for k = 1:size(degree,2)
                    structure_new = rotstr(structure_tran,axis2,degree(k));
                    [protein_starterA,protein_starterB,repulse_num] = pocket2find_PL_Dock(pocket_new,structure_new,RcutoffPL);
                    [Part_Matr_com1_new,CNN1_new,VN1_new] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d); %inter Protein-ligand interactions
                    ZE_part_com1_new = -0.5918*Part_Matr_com1_new - VN1_new*-0.5918;
                    ZN_part_com1_new = ((sqrt(CNN1_new*8+1)+1)/2-4)*log(2)+((sqrt(CNN1_new*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1_new*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1_new*8+1)+1)/2-1)*log(const);
                    ligand_new_all(ln).dG = (ZE_part_com1_new + ZN_part_com1_new*-0.5918)/20;
                    ligand_new_all(ln).coordinate = structure_new;
                    ligand_new_all(ln).repulse_num = repulse_num;
                    ln=ln+1;
                end
            end
        end
    end
    
    for i = 1:size(ligand_new_all,2)
        if ~isempty(ligand_new_all(i).dG)
            dG_select(i,1)=ligand_new_all(i).dG;
            dG_select(i,2)=ligand_new_all(i).repulse_num;
        elseif isempty(ligand_new_all(i).dG)
            dG_select(i,1)=100;
            dG_select(i,2)=ligand_new_all(i).repulse_num;
        end
    end
    
    dG_select(:,3)=1:size(dG_select,1);
    
    dG_select(dG_select(:,2)>=10,:)=[];
    if size(dG_select,1)==0
        break
    end
    dG_select = sortrows(dG_select,[1]);
    
    dG_select_pass=dG_select(1,:);
    
    
    if dG_select_pass(1,1)-E_evaluate(1,1)<=-0.01
        
        ligand_input = ligand_new_all(dG_select_pass(1,3)).coordinate;
        repulse_num_loop = ligand_new_all(dG_select_pass(1,3)).repulse_num;
        loop=loop+1;
        
    elseif (dG_select_pass(1,1)-E_evaluate(1,1)>-0.01 | repulse_num_loop == 0)
        break
    end
    
end
ligand_new_str_structure = ligand_input;
ligand_new_str_dG = dG_select_pass;
ligand_new_str_repulse_num = repulse_num_loop;


