function com_vdw_new = BF_convert(com_vdw_new, ligand, jz)

ring_mark = ligand(jz,8);
ring_size = ligand(jz,9);
C2ar_c3adj=[];
Car_c3adj=[];
C3_atom=[];
C2_atom=[];
C2X_atom=[];
CarX_atom=[];
C2arX_atom=[];
Car_candi=[];
C2ar_candi=[];
O3_atom=[];
O2_atom=[];
O3H_atom=[];
Oall_atom=[];
NH_atom=[];
Nall_atom=[];
Nar_atom=[];
C2Car_atom=[];
ON_atom = [];
Polar_Aro_atom = [];

R1=1.9;
R2=2.1;
R3=4.5;
R4=6;
ad=1;

C3_atom = ligand((ligand(:,4) == 1 | ligand(:,4) == 27) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);

C2X_atom = ligand(ligand(:,4) == 30 & dist(ligand,ligand(jz,:))<=R4 & dist(ligand,ligand(jz,:))>1,:);
CarX_atom = ligand(ligand(:,4) == 33 & dist(ligand,ligand(jz,:))<=R4 & dist(ligand,ligand(jz,:))>1,:);
C2arX_atom = ligand((ligand(:,4) == 30|ligand(:,4) == 33) & dist(ligand,ligand(jz,:))<=R4 & dist(ligand,ligand(jz,:))>1,:);

Car_candi = ligand(ligand(:,4) == 4 | ligand(:,4) == 33,:);
C2ar_candi = ligand(ligand(:,4) == 2 | ligand(:,4) == 4 | ligand(:,4) == 30 | ligand(:,4) == 33,:);

O3_atom = ligand((ligand(:,4) == 5 | ligand(:,4) == 6) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
O2_atom = ligand((ligand(:,4) >= 7 & ligand(:,4) <= 9) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Oall_atom = ligand((ligand(:,4) >= 5 & ligand(:,4) <= 9) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);


NH_atom = ligand((ligand(:,4) >= 11 & ligand(:,4) <= 13) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Nall_atom = ligand((ligand(:,4) >= 10 & ligand(:,4) <= 13) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Nar_atom = ligand(ligand(:,4) == 10 & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
C2Car_atom = ligand((ligand(:,4) == 2 | ligand(:,4) == 4 | ligand(:,4) == 30 | ligand(:,4) == 33) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Call_atom = ligand((ligand(:,4) == 1 |ligand(:,4) == 2 | ligand(:,4) == 4 | ligand(:,4) == 27 | ligand(:,4) == 30 | ligand(:,4) == 33) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Halo_atom = ligand((ligand(:,4) == 15 |ligand(:,4) == 16 | ligand(:,4) == 17 | ligand(:,4) == 18) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
S_atom = ligand((ligand(:,4) == 20 | ligand(:,4) == 21) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
SP_atom = ligand((ligand(:,4) == 20 | ligand(:,4) == 21 | ligand(:,4) == 14) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
ON_atom = ligand((ligand(:,4) >= 5 & ligand(:,4) <= 13) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Polar_Aro_atom = ligand(((ligand(:,4) >= 5 & ligand(:,4) <= 13)|(ligand(:,4) >= 15 & ligand(:,4) <= 18)) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);

H_donor_atom = ligand((ligand(:,4) == 5 | ligand(:,4) == 11 | ligand(:,4) == 12 | ligand(:,4) == 13) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
H_acceptor_atom = ligand((ligand(:,4) == 5 | ligand(:,4) == 7 | ligand(:,4) == 8 | ligand(:,4) == 9) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
C2ar_atom = ligand((ligand(:,4) == 2 | ligand(:,4) == 4 | ligand(:,4) == 30 | ligand(:,4) == 33) & ligand(:,7) ~= ligand(jz,7),:);
Car_atom = ligand((ligand(:,4) == 4 | ligand(:,4) == 33) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
Other_atom = ligand(ligand(:,7) ~= ligand(jz,7) & dist(ligand,ligand(jz,:))<=R3,:);


for i = 1:size(C2ar_candi,1)
    connect_num = connectdetect2(C2ar_candi(i,:),C3_atom,R1);
    if connect_num==1
        ligand(ligand(:,7)==C2ar_candi(i,7),12)=1;
        C2ar_c3adj(ad,:)=C2ar_candi(i,:);
        ad=ad+1;
    elseif connect_num>1
        ligand(ligand(:,7)==C2ar_candi(i,7),12)=2;
    else
        ligand(ligand(:,7)==C2ar_candi(i,7),12)=0;
    end
end


ligand(ligand(:,4)==30,12)=0;










if ligand(jz,4) == 1 || ligand(jz,4) == 27
    %skhs1 = skhs1^(1/-180);
    com_vdw_new = com_vdw_new.^(1/-180);
elseif ligand(jz,4) == 2
    %skhs1 = skhs1^(1/-132);
    com_vdw_new = com_vdw_new.^(1/-132);
    
    
    connect_num = connectdetect(C2_atom);
    if connect_num>=2 && vectofmat(ligand(jz,:),C2_atom)==1
        %skhs1 = skhs1^(1/3.6);
        com_vdw_new = com_vdw_new.^(1/3.6);
    end
    
    
    C3_atom = ligand(ligand(:,4) == 1,:);
    connect_num = connectdetect2(ligand(jz,:),C3_atom,R1);
    if connect_num>=1
        %skhs1 = skhs1^(1/4);
        com_vdw_new = com_vdw_new.^(1/4);
    end
    
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(-1/(1.803*0.8));
        com_vdw_new = com_vdw_new.^(-1/(1.803*0.8));
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(-1/(1.803*0.8));
        com_vdw_new = com_vdw_new.^(-1/(1.803*0.8));
    end
    
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,SP_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(-1/(1.803*0.9));
        com_vdw_new = com_vdw_new.^(-1/(1.803*0.9));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
    if torsion_mark>=1
        %skhs1 = skhs1^(1/-1.803);
        com_vdw_new = com_vdw_new.^(1/-1.803);
    end
    
    
    torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2arX_atom);
    if torsion_mark>=1
        %skhs1 = skhs1^(1/(-1.803));
        com_vdw_new = com_vdw_new.^(1/(-1.803));
    end
    
elseif ligand(jz,4) ==3
    %skhs1 = skhs1^(1/106);
    com_vdw_new = com_vdw_new.^(1/106);
    connect_num = connectdetect2(ligand(jz,:),Other_atom,R1);
    if connect_num == 1
        %skhs1 = skhs1^(1/3);
        com_vdw_new = com_vdw_new.^(1/3);
    end
    connect_numN = connectdetect2(ligand(jz,:),Nall_atom,R1);
    if connect_num == 2 && connect_numN==1
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
    
elseif ligand(jz,4) ==30
    %skhs1 = skhs1^(1/238);
    com_vdw_new = com_vdw_new.^(1/238);
    ligand_other = ligand(ligand(:,7) ~= ligand(jz,7),:);
    connect_num1 = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    connect_num2 = connectdetect2(ligand(jz,:),ligand_other,R1);
    if (connect_num1>=1 || connect_num2 >= 3)
        %skhs1 = skhs1^(1/(0.6*0.2));
        com_vdw_new = com_vdw_new.^(1/(0.6*0.2));
    elseif (connect_num1==0 || connect_num2 == 2)
        %skhs1 = skhs1^(1/(0.23));
        com_vdw_new = com_vdw_new.^(1/(0.23));
    end
    connect_num = connectdetect2(ligand(jz,:),Oall_atom,R1);
    if connect_num==2
        %skhs1 = skhs1^(1/2.3);
        com_vdw_new = com_vdw_new.^(1/2.3);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%angle_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
    %%%if angle_mark>=1
    %%%    %skhs1 = skhs1^(1/(0.8));
    %%%    com_vdw_new = com_vdw_new.^(1/(0.8));
    %%%end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
    if torsion_mark>=1
        %skhs1 = skhs1^(1/(0.8));
        com_vdw_new = com_vdw_new.^(1/(0.8));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/(0.8*1));
        com_vdw_new = com_vdw_new.^(1/(0.8*1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.8*1));
        com_vdw_new = com_vdw_new.^(1/(0.8*1));
    end
    
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,SP_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.9*1));
        com_vdw_new = com_vdw_new.^(1/(0.9*1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2arX_atom);
    if torsion_mark>=1
        %skhs1 = skhs1^(1/(0.8*1));
        com_vdw_new = com_vdw_new.^(1/(0.8*1));
    end
    
    
elseif ligand(jz,4) == 4 || ligand(jz,4) == 33
    %skhs1 = skhs1^(1/238);
    com_vdw_new = com_vdw_new.^(1/238);
    
    %Car_atom = ligand(ligand(:,4) == 4 & ligand(:,7) ~= ligand(jz,7),:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ligand(jz,12) ==1
        %skhs1 = skhs1^(1/0.7);
        com_vdw_new = com_vdw_new.^(1/0.7);
    end
    
    %connect_numOH = connectdetect2(ligand(jz,:),OH_atom);
    
    
    %ligand_other = ligand(ligand(:,7) ~= ligand(jz,7),:);
    %connect_num1 = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    %connect_num2 = connectdetect2(ligand(jz,:),ligand_other,R1);
    connect_num = connectdetect2(ligand(jz,:),Car_atom,R1);
    if connect_num==3
        %skhs1 = skhs1^(1/0.2);
        com_vdw_new = com_vdw_new.^(1/0.2);
    end
    
    
    %if (ligand(jz,4) ==33|ligand(jz,4) ==30) && (connect_num1>=1 | connect_num2 >= 3)
    if (ligand(jz,4) ==33)
        %skhs1 = skhs1^(1/(0.7*1));
        com_vdw_new = com_vdw_new.^(1/(0.7*1));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,C3_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.8));
        com_vdw_new = com_vdw_new.^(1/0.8);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    angle_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.8));
        com_vdw_new = com_vdw_new.^(1/(0.8));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
    if torsion_mark>=1
        %skhs1 = skhs1^(1/(0.8));
        com_vdw_new = com_vdw_new.^(1/(0.8));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%angle_mark = angle_find_W(ligand(jz,:),Car_atom,O3_atom);
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.8*1));
        com_vdw_new = com_vdw_new.^(1/(0.8*1));
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.8*1));
        com_vdw_new = com_vdw_new.^(1/(0.8*1));
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,SP_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/(0.9*1));
        com_vdw_new = com_vdw_new.^(1/(0.9*1));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2arX_atom);
    if torsion_mark>=1
        %skhs1 = skhs1^(1/(0.8*1));
        com_vdw_new = com_vdw_new.^(1/(0.8*1));
    end
    
    
elseif ligand(jz,4) == 10
    %skhs1 = skhs1^(1/22);
    com_vdw_new = com_vdw_new.^(1/22);
    connect_num = connectdetect2(ligand(jz,:),Nall_atom,R1);
    if connect_num>=1
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    %%%%%%%%%%%%%%%%%%%judge Nar in an aromatic ring%%%%%%%%%%%
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    if connect_num==2
        %skhs1 = skhs1^(1/0.8);
        com_vdw_new = com_vdw_new.^(1/0.8);
        para_atom_final=[];
        para_atom_final = para_atom_find(ligand(jz,:),C2ar_atom,Nar_atom);
        
        
        if ~isempty(para_atom_final)
            %skhs1 = skhs1^(1/(1.6));
            com_vdw_new = com_vdw_new.^(1/(1.6));
            %break
        end
        
        
        angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
        if angle_mark>=2
            %skhs1 = skhs1^(1/(2));
            com_vdw_new = com_vdw_new.^(1/(2));
        end
        
        
    end
    
    
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/1.6);
        com_vdw_new = com_vdw_new.^(1/1.6);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/1.6);
        com_vdw_new = com_vdw_new.^(1/1.6);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,C2ar_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/1);
        com_vdw_new = com_vdw_new.^(1/1);
    end
    
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
elseif ligand(jz,4) >= 11 && ligand(jz,4) <= 13
    %skhs1 = skhs1^(1/19);
    com_vdw_new = com_vdw_new.^(1/19);
    
    %N3H_atom = ligand(ligand(:,4) >= 11 & ligand(:,4) <= 13,:);
    connect_num = connectdetect2(ligand(jz,:),Nall_atom,R1);
    if connect_num>=1
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
    if angle_mark>=2
        
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    
    C2arXjz_atom = [];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
    H_acceptor_atom_input = [];
    H_acceptor_atom_input = H_acceptor_atom(dist(H_acceptor_atom,ligand(jz,:))<3.0,:);
    intra_H_bond_atom = [];
    intra_H_bond_atom = H_acceptor_atom_input(dist(H_acceptor_atom_input,ligand(jz,:))<3.0 & dist(H_acceptor_atom_input,ligand(jz,:))>2.55,:);
    
    if ~isempty(intra_H_bond_atom)
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    
    
    
    
    
elseif ligand(jz,4) == 5
    %skhs1 = skhs1^(1/19);
    com_vdw_new = com_vdw_new.^(1/19);
    connect_num = connectdetect2(ligand(jz,:),O3_atom,R1);
    if connect_num>=1
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    if connect_num>=1
        
        %skhs1 = skhs1^(1/1);
        com_vdw_new = com_vdw_new.^(1/1);
    end
    
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    if connect_num>=1
        
        %skhs1 = skhs1^(1/1);
        com_vdw_new = com_vdw_new.^(1/1);
    end
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,O2_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/0.9);
        com_vdw_new = com_vdw_new.^(1/0.9);
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/0.9);
        com_vdw_new = com_vdw_new.^(1/0.9);
    end
    
    
    %connect_num = connectdetect2(ligand(jz,:),C2arX_atom,R1);
    %if connect_num>=1
    %    %skhs1 = skhs1^(1/0.9);
    %    com_vdw_new = com_vdw_new.^(1/0.9);
    %end
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
        
    end
    H_acceptor_atom_input = [];
    H_acceptor_atom_input = H_acceptor_atom(dist(H_acceptor_atom,ligand(jz,:))<3.0,:);
    intra_H_bond_atom1 = [];
    intra_H_bond_atom1 = H_acceptor_atom_input(dist(H_acceptor_atom_input,ligand(jz,:))<3.0 & dist(H_acceptor_atom_input,ligand(jz,:))>2.55,:);
    
    H_donor_atom_input = [];
    H_donor_atom_input = H_donor_atom(dist(H_donor_atom,ligand(jz,:))<3.0,:);
    intra_H_bond_atom2 = [];
    intra_H_bond_atom2 = H_donor_atom_input(dist(H_donor_atom_input,ligand(jz,:))<3.0 & dist(H_donor_atom_input,ligand(jz,:))>2.55,:);
    
    if ~isempty(intra_H_bond_atom1) || ~isempty(intra_H_bond_atom2)
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    
elseif ligand(jz,4) == 6
    %skhs1 = skhs1^(1/35);
    com_vdw_new = com_vdw_new.^(1/35);
    connect_num = connectdetect2(ligand(jz,:),O3_atom,R1);
    if connect_num>=1
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
    if connect_num>=1
        
        %skhs1 = skhs1^(1/1);
        com_vdw_new = com_vdw_new.^(1/1);
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,O2_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/1.4);
        com_vdw_new = com_vdw_new.^(1/1.4);
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/1.4);
        com_vdw_new = com_vdw_new.^(1/1.4);
    end
    
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
elseif ligand(jz,4) >= 7 && ligand(jz,4) <= 9
    %skhs1 = skhs1^(1/28);
    com_vdw_new = com_vdw_new.^(1/28);
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/1.4);
        com_vdw_new = com_vdw_new.^(1/1.4);
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
    if angle_mark>=1
        %skhs1 = skhs1^(1/1.4);
        com_vdw_new = com_vdw_new.^(1/1.4);
    end
    
    angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,C2ar_atom);
    if angle_mark>=1
        
        %skhs1 = skhs1^(1/2.1);
        com_vdw_new = com_vdw_new.^(1/2.1);
    end
    
    if ligand(jz,4) == 9
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    connect_num1 = connectdetect2(ligand(jz,:),Nar_atom,R1);
    angle_mark1 = angle_find_W(ligand(jz,:),Nar_atom,O2_atom);
    if angle_mark1>=1 && connect_num1>=1
        
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
    
    
    
    
    
    H_donor_atom_input = [];
    H_donor_atom_input = H_donor_atom(dist(H_donor_atom,ligand(jz,:))<3.0,:);
    intra_H_bond_atom = [];
    intra_H_bond_atom = H_donor_atom_input(dist(H_donor_atom_input,ligand(jz,:))<3.0 & dist(H_donor_atom_input,ligand(jz,:))>2.55,:);
    
    if ~isempty(intra_H_bond_atom)
        %skhs1 = skhs1^(1/2);
        com_vdw_new = com_vdw_new.^(1/2);
    end
    %%%torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_atom);
    %%%if torsion_mark>=1
    %%%    %skhs1 = skhs1^(1/tr2);
    %%%    com_vdw_new = com_vdw_new.^(1/tr2);
    %%%end
    
elseif ligand(jz,4) == 20 || ligand(jz,4) == 21 || ligand(jz,4) == 14
    
    
    %skhs1 = skhs1^(1/36);
    com_vdw_new = com_vdw_new.^(1/36);
    
    connect_num = connectdetect2(ligand(jz,:),Other_atom,R2);
    connect_numS = connectdetect2(ligand(jz,:),SP_atom,R2);
    if connect_numS >= 1
        %skhs1 = skhs1^(1/1.4);
        com_vdw_new = com_vdw_new.^(1/1.4);
    elseif connect_numS == 0 && connect_num >= 2
        %skhs1 = skhs1^(1/1);
        com_vdw_new = com_vdw_new.^(1/1);
    end
    
    
    
    
    
    
    
    
    
    
    %angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
    
    
elseif ligand(jz,4) == 17 || ligand(jz,4) == 18
    %skhs1 = skhs1^(1/67);
    com_vdw_new = com_vdw_new.^(1/67);
    angle_mark = angle_find_W(ligand(jz,:),C3_atom,Halo_atom);
    if angle_mark>=2
        %skhs1 = skhs1^(1/3.7);
        com_vdw_new = com_vdw_new.^(1/3.7);
    end
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R2);
    if connect_num>=1
        %skhs1 = skhs1^(1/3.7);
        com_vdw_new = com_vdw_new.^(1/3.7);
    end
    
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
    
elseif ligand(jz,4) == 15
    %skhs1 = skhs1^(-1/81);
    com_vdw_new = com_vdw_new.^(-1/81);
    angle_mark = angle_find_W(ligand(jz,:),C3_atom,Halo_atom);
    if angle_mark>=2
        %skhs1 = skhs1^(1/3.7);
        com_vdw_new = com_vdw_new.^(1/3.7);
    end
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R2);
    if connect_num>=1
        %skhs1 = skhs1^(1/3.7);
        com_vdw_new = com_vdw_new.^(1/3.7);
    end
    
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
elseif ligand(jz,4) == 16
    %skhs1 = skhs1^(1/119);
    com_vdw_new = com_vdw_new.^(1/119);
    angle_mark = angle_find_W(ligand(jz,:),C3_atom,Halo_atom);
    if angle_mark>=2
        %skhs1 = skhs1^(1/3.7);
        com_vdw_new = com_vdw_new.^(1/3.7);
    end
    connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R2);
    if connect_num>=1
        %skhs1 = skhs1^(1/3.7);
        com_vdw_new = com_vdw_new.^(1/3.7);
    end
    C2arXjz_atom=[];
    oppo_atom_final1 = [];
    oppo_atom_final2 = [];
    C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
    [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
    if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
        %skhs1 = skhs1^(1/1.5);
        com_vdw_new = com_vdw_new.^(1/1.5);
    end
end

if ligand(jz,4) == 27
    %skhs1 = skhs1^(1/-1.6);
    com_vdw_new = com_vdw_new.^(1/-1.6);
end
%if ligand(jz,4) == 30
%    %skhs1 = skhs1^(1/0.2);
%    com_vdw_new = com_vdw_new.^(1/0.2);
%end

%%%cf_vdw_new=sum(cf_vdw);
if ring_mark == 1 && ring_size == 3 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
    %if ring_mark == 1 && ring_size == 3 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
    %skhs1 = skhs1^(1/2.8);
    com_vdw_new = com_vdw_new^(1/2.8);
    %%skhs1 = skhs1^(1/3);
    %com_vdw_new = com_vdw_new^(1/3);
elseif ring_mark == 1 && ring_size == 4 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
    %elseif ring_mark == 1 && ring_size == 4 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
    %skhs1 = skhs1^(1/2.4);
    com_vdw_new = com_vdw_new^(1/2.4);
    %%skhs1 = skhs1^(1/2.8);
    %com_vdw_new = com_vdw_new^(1/2.8);
elseif ring_mark == 1 && ring_size == 5 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
    %elseif ring_mark == 1 && ring_size == 5 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
    %skhs1 = skhs1^(1/2.2);
    com_vdw_new = com_vdw_new^(1/2.2);
    %%skhs1 = skhs1^(1/2.6);
    %com_vdw_new = com_vdw_new^(1/2.6);
elseif ring_mark == 1 && ring_size == 6 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
    %elseif ring_mark == 1 && ring_size == 6 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
    %skhs1 = skhs1^(1/2);
    com_vdw_new = com_vdw_new^(1/2);
    %%skhs1 = skhs1^(1/2);
    %com_vdw_new = com_vdw_new^(1/2);
end