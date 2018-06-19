clear all
home_dir = strcat('input structure file directory\');
out_dir = strcat('MTDocking output folder\');
list1=dir(home_dir);
size_list1=size(list1);
dist_cent_search = 4;
RcutoffPL=6;
R_srch=0.5;
MTsample = R_srch;
R_L_tor_srch=0.1;
rd = 0.005;
hnstate = R_srch/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2;


Rpocketcutoff=5;
tn=0;
clear test_name
code_dir = pwd;
jpz=1;


for zh=1:size_list1(1,1)
    tic
    if length(list1(zh,1).name)<8
        continue
    elseif length(list1(zh,1).name)>=8
        list1_name=list1(zh,1).name;
       
        if strcmp(list1_name(end-3:end),'.pdb')==1
            protein_id = list1_name(1:end-4);
            pdbid = protein_id;
            
            protein_name = list1_name;
            protein_dir = [ strcat(home_dir,protein_name)];
            [protein,protein_Namelist] = protein_define_final(protein_dir);
            
            
            for hl=1:size_list1(1,1)
                
                if length(list1(hl,1).name)<8
                    continue
                elseif length(list1(hl,1).name)>=8
                    
                    
                    list1_name1=list1(hl,1).name;
                    if strcmp(list1_name1(end-4:end),'.mol2')==1
                        
                        ligand_name = list1_name1;
                        ligand_dir = [ strcat(home_dir,ligand_name)];
                        
                        
                        
                        [ligand,name_str]= ligand_define_sol_spec(ligand_dir);
                        
                        [pocket] = pocket2find_PL_region(protein,ligand,Rpocketcutoff);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        output_dir = strcat(out_dir,protein_name,'\',ligand_name,'\output\');
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % create the heatmap within the binding pocket
                        clear Grid_Heatmap
                        Grid_Heatmap = Heatmap_2a(pocket,protein,ligand);
                        %
                        
                        toc
                        
                        run GARF_Potential
                        
                        tic
                        clear pocket_centroid
                        pocket_centroid(1,1) = mean(pocket(:,1));
                        pocket_centroid(1,2) = mean(pocket(:,2));
                        pocket_centroid(1,3) = mean(pocket(:,3));
                        
                        pocket_centroid = repmat(pocket_centroid,7,1);
                        for i = 1:6
                            dan=ceil(i/2);
                            pocket_centroid(1+i,dan) = pocket_centroid(1,dan)+dist_cent_search*(-1)^i;
                        end
                        
                        
                        
                        
                        ligand_centroid(1,1)=mean(ligand(:,1));
                        ligand_centroid(1,2)=mean(ligand(:,2));
                        ligand_centroid(1,3)=mean(ligand(:,3));
                        
                        ligand_sig = ligandsig(ligand);
                        
                        % place ligand into the mt-scored heatmap grid points
                        [ligand_final] = DioP_anchor_v3(ligand,Grid_Heatmap,pocket_centroid,pocket,protein);
                        %
                        toc
                        
                        clear ligand_final_score
                        ligand_final_score=zeros(size(ligand_final,2),2);
                        for lf=1:size(ligand_final,2)
                            ligand_final_score(lf,1) = ligand_final(lf).dG;
                            ligand_final_score(lf,2) = lf;
                        end
                        
                        ligand_final_score = sortrows(ligand_final_score,[1]);
                        clear ligand_select_X1
                        for lsxn=1:size(ligand_final,2)
                            
                            ligand_select_X1(lsxn).coordinate = ligand_final(ligand_final_score(lsxn,2)).structure;
                            ligand_select_X1(lsxn).dG_X1 = ligand_final(ligand_final_score(lsxn,2)).dG;
                            ligand_select_X1(lsxn).repulse = ligand_final(ligand_final_score(lsxn,2)).repulse;
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        clear dG_select_X1 dG_select_X2 dG_select_pass_X1 dG_select_pass_X2 Final_result_X1 Final_result_X2
                        
                        for i = 1:size(ligand_select_X1,2)
                            if ~isempty(ligand_select_X1(i).dG_X1)
                                if ligand_select_X1(i).repulse<=60
                                    dG_select_X1(i,1)=ligand_select_X1(i).dG_X1;
                                    dG_select_X1(i,2)=ligand_select_X1(i).repulse;
                                elseif ligand_select_X1(i).repulse>60
                                    dG_select_X1(i,1)=0;
                                    dG_select_X1(i,2)=0;
                                end
                            elseif isempty(ligand_select_X1(i).dG_X1)
                                dG_select_X1(i,1)=0;
                                dG_select_X1(i,2)=0;
                            end
                        end
                        
                        dG_select_X1(:,3)=1:size(dG_select_X1,1);
                        dG_select_X1(dG_select_X1(:,1)==0,:)=[];
                        
                        dG_select_X1 = sortrows(dG_select_X1,[1]);
                        
                        if size(dG_select_X1,1)>5
                            dG_select_pass_X1 = dG_select_X1(1:5,:);
                        else
                            dG_select_pass_X1 = dG_select_X1;
                        end
                        
                        dG_select_pass_X1(dG_select_pass_X1(:,1)==0,:)=[];
                        
                        dG_select_pass_X1 = sortrows(dG_select_pass_X1,[1]);
                        
                        % energy relaxation using the function ErelaxB
                        
                        mkdir(output_dir);
                        for i = 1:size(dG_select_pass_X1,1)
                            [ligand_final_str_X1,ligand_final_dG_X1,ligand_final_repulse_num_X1] = ErelaxB(ligand_select_X1,dG_select_pass_X1,i,pocket,vdw_d);
                            
                            RMSD_final = RMSD_cal(ligand_final_str_X1,ligand);
                            
                            Final_result_X1(i,1)=ligand_final_dG_X1(1,1);
                            Final_result_X1(i,2)=ligand_final_repulse_num_X1;
                            Final_result_X1(i,3)=RMSD_final;
                            
                            input_struc = ligand_final_str_X1;
                            listC = strcat(ligand_name,'_pose');
                            
                            out_X1 = lig_output(input_struc,listC,i,output_dir);
                            
                            
                        end
                        
                        
                        %Final_result=[Final_result_X1
                        %    ];
                        
                        %data_output(Final_result,output_dir,code_dir);
                        
                        
                        %Final_result_data(jpz).data=Final_result;
                        %jpz=jpz+1;
                        
                    end
                end
            end
            
        end
    end
end

