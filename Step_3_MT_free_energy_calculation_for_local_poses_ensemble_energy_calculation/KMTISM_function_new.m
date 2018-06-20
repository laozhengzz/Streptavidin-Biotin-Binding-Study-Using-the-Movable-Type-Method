function [Part_Matr_san1,Part_Matr_san2,SA_list_sum] = KMTISM_function_new(protein)
GARF_Potential_OH_ligand_contact
tic

%protein = protein_sol;
%water_d = vdw_d;
%clear vdw_d

rd = 0.005;
halv=100;
halv1=10;

const=halv*2;
const1=halv1*2;
column_num=100;

op_radi=[0
    3.325
    3.73
    3.65
    3.635
    3.66
    3.155
    2.895
    2.85
    2.955
    2.785
    2.84
    2.845
    2.735
    2.82
    2.745
    3.015
    2.76
    3.78
    ];
%list=flipud(list);

C=1  ;
Cs=2 ;
CA=3 ;
CB=4 ;
CC=5 ;
CN=6 ;
CR=7 ;
CT=8 ;
CV=9 ;
CW=10;
H=11 ;
HO=12;
N=13;
N2=14;
N3=15;
NA=16;
NB=17;
O=18 ;
O2=19;
OH=20;
S=21;
SH=22;



C_3=1;
C_2=2;
C_1=3;
C_ar=4;
O_3=5;
O_3p=6;
O_2=7;
O_co2=8;
O_2v=9;
N_2=10;
N_am=11;
N_pl3=12;
N_4=13;
P=14;
F=15;
Cl=16;
Br=17;
I=18;
C_cat=19;
S_3=20;
S_o=21;
HNa=22;
HOa=23;
HSa=24;
H3 =25;
C_3Oa=26;
C_3Na=27;
C_3La=28;
C_2Oa=29;
C_2Na=30;
C_2La=31;
C_arOa=32;
C_arNa=33;
C_arLa=34;
HCa=35;

N_3=12;
N_1=10;
N_ar=10;
S_2=21;
S_o2=22;





K    = 29 ;
MG   = 30  ;
CAA   = 31	;
AL   = 32	;
MN   = 33  	;
FE   = 34    ;
CO   = 35    ;
NI   = 36    ;
CU   = 37    ;
ZN   = 38    ;


N_A=0;
D=1;
D2=2;
A=3;
DA=4;

ALA=101;
ARG=102;
ASN=103;
ASP=104;
CYS=105;
GLN=106;
GLU=107;
GLY=108;
HIS=109;
ILE=110;
LEU=111;
LYS=112;
MET=113;
PHE=114;
PRO=115;
SER=116;
THR=117;
TRP=118;
TYR=119;
VAL=120;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pl_complex=[protential_protein
%    ligand];

%size_pl_complex=size(pl_complex);

%complex_atom_num=[1:1:size_pl_complex(1,1)];

%pl_complex(:,7)=transpose(complex_atom_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ligand=protein;
size_ligand=size(ligand);




ao=-0.3;
bo=0.2;
co=1;
ddo=7.3;
eo=9.8;
fo=6.2;
go=6;
ho=10.1;
jo=-0.5;
ko=0.4;
lo=0.6;
mo=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%water model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Car_candi = ligand(ligand(:,4) == 4 | ligand(:,4) == 33,:);
C2ar_candi = ligand(ligand(:,4) == 2 | ligand(:,4) == 4 | ligand(:,4) == 30 | ligand(:,4) == 33,:);
    

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

%Emult=[];
%28
%0.4
Part_Matr_san1 = 0;
Part_Matr_san2 = 0;
%for tr3=1:0.1:1

ligand_grid=[];
size_ligand_grid=[0,0];
aj=1;
clear vn
vn = zeros(size_ligand(1,1),1);
SA_expo_list = [];
SA_total = 0;
SA_list = [];
Freq_Matr_san=0;



clear grid_coor_stand
grid_d=0.2;
max_x=0+6;
min_x=0-6;
max_y=0+6;
min_y=0-6;
max_z=0+6;
min_z=0-6;
sig_x=transpose(min_x:grid_d:max_x);
sig_y=transpose(min_y:grid_d:max_y);
sig_z=transpose(min_z:grid_d:max_z);

size_sig_x=size(sig_x);
size_sig_y=size(sig_y);
size_sig_z=size(sig_z);


for q=1:1:size_sig_x(1,1)
    grid_coor_stand((q-1)*size_sig_y(1,1)*size_sig_z(1,1)+1:q*size_sig_y(1,1)*size_sig_z(1,1),1)=sig_x(q,1);
end
%toc
%tic
for q=1:1:size_sig_y(1,1)
    grid_coor1r((q-1)*size_sig_z(1,1)+1:q*size_sig_z(1,1),1)=sig_y(q,1);
end
%toc
%tic
grid_coor_stand(:,2)=repmat(grid_coor1r,size_sig_x(1,1),1);
grid_coor_stand(:,3)=repmat(sig_z,size_sig_x(1,1)*size_sig_y(1,1),1);

toc



Part_Matr_san=0;
for jz=1:1:size_ligand(1,1)

    ring_mark = ligand(jz,8);
    ring_size = ligand(jz,9);
    
    C3_atom = ligand((ligand(:,4) == 1 | ligand(:,4) == 27) & dist(ligand,ligand(jz,:))<=R3 & dist(ligand,ligand(jz,:))>1,:);
    
    C2X_atom = ligand(ligand(:,4) == 30 & dist(ligand,ligand(jz,:))<=R4 & dist(ligand,ligand(jz,:))>1,:);
    CarX_atom = ligand(ligand(:,4) == 33 & dist(ligand,ligand(jz,:))<=R4 & dist(ligand,ligand(jz,:))>1,:);
    C2arX_atom = ligand((ligand(:,4) == 30|ligand(:,4) == 33) & dist(ligand,ligand(jz,:))<=R4 & dist(ligand,ligand(jz,:))>1,:);
    
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
    
    
    
    
    clear grid_coor
    grid_coor(:,1) = grid_coor_stand(:,1) + ligand(jz,1);
    grid_coor(:,2) = grid_coor_stand(:,2) + ligand(jz,2);
    grid_coor(:,3) = grid_coor_stand(:,3) + ligand(jz,3);
    %C2ar_atom = ligand((ligand(:,4) == 2 | ligand(:,4) == 4 | ligand(:,4) == 30 | ligand(:,4) == 33) & ligand(:,7) ~= ligand(jz,7),:);
    %Car_atom = ligand((ligand(:,4) == 4 | ligand(:,4) == 33) & ligand(:,7) ~= ligand(jz,7),:);
    %Other_atom = ligand(ligand(:,7) ~= ligand(jz,7),:);
    for lyr = ligand(jz,6):grid_d:6+grid_d
        clear P_sh
        P_sh=grid_coor(  (sqrt( (grid_coor(:,1)-ligand(jz,1)).^2+(grid_coor(:,2)-ligand(jz,2)).^2+(grid_coor(:,3)-ligand(jz,3)).^2 )<= lyr+grid_d)&(sqrt( (grid_coor(:,1)-ligand(jz,1)).^2+(grid_coor(:,2)-ligand(jz,2)).^2+(grid_coor(:,3)-ligand(jz,3)).^2 )> lyr),:);
        
        R_point=ligand(jz,:);
        lig_blackwhole=ligand(sqrt( (ligand(:,1)-R_point(1,1)).^2+(ligand(:,2)-R_point(1,2)).^2+(ligand(:,3)-R_point(1,3)).^2 )<7.8,:);
        
        size_lig_blackwhole=size(lig_blackwhole);
        
        
        
        for jj=1:size_lig_blackwhole(1,1)
            if lig_blackwhole(jj,7)~=ligand(jz,7)
                P_sh(sqrt( (P_sh(:,1)-lig_blackwhole(jj,1)).^2+(P_sh(:,2)-lig_blackwhole(jj,2)).^2+(P_sh(:,3)-lig_blackwhole(jj,3)).^2 )<=lig_blackwhole(jj,6),:)=[];
            end
        end
        size_P_sh=size(P_sh);
        
        SA_cho=size_P_sh(1,1)*grid_d*grid_d;
        
        %SA_list(jz,1) = SA_cho;
        
        if lyr == ligand(jz,6)
            SA_list(jz,1) = SA_cho;
        end
        skhs1=exp(SA_cho*0.034165);
        skhs2=skhs1;
        %vn(jz,1) = log(skhs1);
        sd=find(abs(water_d(:,1)-ligand(jz,6))<=rd/2);
        com_vdw=water_d(sd-rd/0.01:sd+rd/0.01,ligand(jz,4)+1);
        
        
        %com_vdw=water_d(sd:1000,ligand(jz,4)+1);
        %%%cf_vdw=water_f(sd:800,ligand(jz,4));
        %%%skhs1=skhs1.^ao;
        
        %%%com_vdw=com_vdw.*cf_vdw;
        
        com_vdw_new=sum(com_vdw);
        com_vdw_new2=com_vdw_new;
        if ligand(jz,4) == 1 || ligand(jz,4) == 27
            skhs1 = skhs1^(1/-180);
            com_vdw_new = com_vdw_new.^(1/-180);
        elseif ligand(jz,4) == 2
            skhs1 = skhs1^(1/-132);
            com_vdw_new = com_vdw_new.^(1/-132);
            
            
            connect_num = connectdetect(C2_atom);
            if connect_num>=2 && vectofmat(ligand(jz,:),C2_atom)==1
                skhs1 = skhs1^(1/3.6);
                com_vdw_new = com_vdw_new.^(1/3.6);
            end
            
            
            C3_atom = ligand(ligand(:,4) == 1,:);
            connect_num = connectdetect2(ligand(jz,:),C3_atom,R1);
            if connect_num>=1
                skhs1 = skhs1^(1/4);
                com_vdw_new = com_vdw_new.^(1/4);
            end
            
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(-1/(1.803*0.8));
                com_vdw_new = com_vdw_new.^(-1/(1.803*0.8));
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(-1/(1.803*0.8));
                com_vdw_new = com_vdw_new.^(-1/(1.803*0.8));
            end
            
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,SP_atom);
            if angle_mark>=1
                skhs1 = skhs1^(-1/(1.803*0.9));
                com_vdw_new = com_vdw_new.^(-1/(1.803*0.9));
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
            if torsion_mark>=1
                skhs1 = skhs1^(1/-1.803);
                com_vdw_new = com_vdw_new.^(1/-1.803);
            end
            
            
            torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2arX_atom);
            if torsion_mark>=1
                skhs1 = skhs1^(1/(-1.803));
                com_vdw_new = com_vdw_new.^(1/(-1.803));
            end
            
        elseif ligand(jz,4) ==3
            skhs1 = skhs1^(1/106);
            com_vdw_new = com_vdw_new.^(1/106);
            connect_num = connectdetect2(ligand(jz,:),Other_atom,R1);
            if connect_num == 1
                skhs1 = skhs1^(1/3);
                com_vdw_new = com_vdw_new.^(1/3);
            end
            connect_numN = connectdetect2(ligand(jz,:),Nall_atom,R1);
            if connect_num == 2 && connect_numN==1
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
            
        elseif ligand(jz,4) ==30
            skhs1 = skhs1^(1/238);
            com_vdw_new = com_vdw_new.^(1/238);
            ligand_other = ligand(ligand(:,7) ~= ligand(jz,7),:);
            connect_num1 = connectdetect2(ligand(jz,:),C2ar_atom,R1);
            connect_num2 = connectdetect2(ligand(jz,:),ligand_other,R1);
            if (connect_num1>=1 || connect_num2 >= 3)
                skhs1 = skhs1^(1/(0.6*0.2));
                com_vdw_new = com_vdw_new.^(1/(0.6*0.2));
            elseif (connect_num1==0 || connect_num2 == 2)
                skhs1 = skhs1^(1/(0.23));
                com_vdw_new = com_vdw_new.^(1/(0.23));
            end
            connect_num = connectdetect2(ligand(jz,:),Oall_atom,R1);
            if connect_num==2
                skhs1 = skhs1^(1/2.3);
                com_vdw_new = com_vdw_new.^(1/2.3);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%angle_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
            %%%if angle_mark>=1
            %%%    skhs1 = skhs1^(1/(0.8));
            %%%    com_vdw_new = com_vdw_new.^(1/(0.8));
            %%%end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
            if torsion_mark>=1
                skhs1 = skhs1^(1/(0.8));
                com_vdw_new = com_vdw_new.^(1/(0.8));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/(0.8*1));
                com_vdw_new = com_vdw_new.^(1/(0.8*1));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.8*1));
                com_vdw_new = com_vdw_new.^(1/(0.8*1));
            end
            
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,SP_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.9*1));
                com_vdw_new = com_vdw_new.^(1/(0.9*1));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2arX_atom);
            if torsion_mark>=1
                skhs1 = skhs1^(1/(0.8*1));
                com_vdw_new = com_vdw_new.^(1/(0.8*1));
            end
            
            
        elseif ligand(jz,4) == 4 || ligand(jz,4) == 33
            skhs1 = skhs1^(1/238);
            com_vdw_new = com_vdw_new.^(1/238);
            
            %Car_atom = ligand(ligand(:,4) == 4 & ligand(:,7) ~= ligand(jz,7),:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ligand(jz,12) ==1
                skhs1 = skhs1^(1/0.7);
                com_vdw_new = com_vdw_new.^(1/0.7);
            end
            
            %connect_numOH = connectdetect2(ligand(jz,:),OH_atom);
            
            
            %ligand_other = ligand(ligand(:,7) ~= ligand(jz,7),:);
            %connect_num1 = connectdetect2(ligand(jz,:),C2ar_atom,R1);
            %connect_num2 = connectdetect2(ligand(jz,:),ligand_other,R1);
            connect_num = connectdetect2(ligand(jz,:),Car_atom,R1);
            if connect_num==3
                skhs1 = skhs1^(1/0.2);
                com_vdw_new = com_vdw_new.^(1/0.2);
            end
            
            
            %if (ligand(jz,4) ==33|ligand(jz,4) ==30) && (connect_num1>=1 | connect_num2 >= 3)
            if (ligand(jz,4) ==33)
                skhs1 = skhs1^(1/(0.7*1));
                com_vdw_new = com_vdw_new.^(1/(0.7*1));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,C3_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.8));
                com_vdw_new = com_vdw_new.^(1/0.8);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            angle_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.8));
                com_vdw_new = com_vdw_new.^(1/(0.8));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_c3adj);
            if torsion_mark>=1
                skhs1 = skhs1^(1/(0.8));
                com_vdw_new = com_vdw_new.^(1/(0.8));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%angle_mark = angle_find_W(ligand(jz,:),Car_atom,O3_atom);
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.8*1));
                com_vdw_new = com_vdw_new.^(1/(0.8*1));
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.8*1));
                com_vdw_new = com_vdw_new.^(1/(0.8*1));
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,SP_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/(0.9*1));
                com_vdw_new = com_vdw_new.^(1/(0.9*1));
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2arX_atom);
            if torsion_mark>=1
                skhs1 = skhs1^(1/(0.8*1));
                com_vdw_new = com_vdw_new.^(1/(0.8*1));
            end
            
            
        elseif ligand(jz,4) == 10
            skhs1 = skhs1^(1/22);
            com_vdw_new = com_vdw_new.^(1/22);
            connect_num = connectdetect2(ligand(jz,:),Nall_atom,R1);
            if connect_num>=1
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            %%%%%%%%%%%%%%%%%%%judge Nar in an aromatic ring%%%%%%%%%%%
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
            if connect_num==2
                skhs1 = skhs1^(1/0.8);
                com_vdw_new = com_vdw_new.^(1/0.8);
                para_atom_final=[];
                para_atom_final = para_atom_find(ligand(jz,:),C2ar_atom,Nar_atom);
                
                
                if ~isempty(para_atom_final)
                    skhs1 = skhs1^(1/(1.6));
                    com_vdw_new = com_vdw_new.^(1/(1.6));
                    break
                end
                
                
                angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
                if angle_mark>=2
                    skhs1 = skhs1^(1/(2));
                    com_vdw_new = com_vdw_new.^(1/(2));
                end
                
                
            end
            
            
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/1.6);
                com_vdw_new = com_vdw_new.^(1/1.6);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/1.6);
                com_vdw_new = com_vdw_new.^(1/1.6);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,C2ar_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/1);
                com_vdw_new = com_vdw_new.^(1/1);
            end
            
            C2arXjz_atom=[];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
        elseif ligand(jz,4) >= 11 && ligand(jz,4) <= 13
            skhs1 = skhs1^(1/19);
            com_vdw_new = com_vdw_new.^(1/19);
            
            %N3H_atom = ligand(ligand(:,4) >= 11 & ligand(:,4) <= 13,:);
            connect_num = connectdetect2(ligand(jz,:),Nall_atom,R1);
            if connect_num>=1
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
            if angle_mark>=2
                
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            
            C2arXjz_atom = [];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
            H_acceptor_atom_input = [];
            H_acceptor_atom_input = H_acceptor_atom(dist(H_acceptor_atom,ligand(jz,:))<3.0,:);
            intra_H_bond_atom = [];
            intra_H_bond_atom = H_acceptor_atom_input(dist(H_acceptor_atom_input,ligand(jz,:))<3.0 & dist(H_acceptor_atom_input,ligand(jz,:))>2.55,:);
            
            if ~isempty(intra_H_bond_atom)
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            
            
            
            
            
        elseif ligand(jz,4) == 5
            skhs1 = skhs1^(1/19);
            com_vdw_new = com_vdw_new.^(1/19);
            connect_num = connectdetect2(ligand(jz,:),O3_atom,R1);
            if connect_num>=1
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
            if connect_num>=1
                
                skhs1 = skhs1^(1/1);
                com_vdw_new = com_vdw_new.^(1/1);
            end
            
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
            if connect_num>=1
                
                skhs1 = skhs1^(1/1);
                com_vdw_new = com_vdw_new.^(1/1);
            end
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,O2_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/0.9);
                com_vdw_new = com_vdw_new.^(1/0.9);
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nar_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/0.9);
                com_vdw_new = com_vdw_new.^(1/0.9);
            end
            
            
            %connect_num = connectdetect2(ligand(jz,:),C2arX_atom,R1);
            %if connect_num>=1
            %    skhs1 = skhs1^(1/0.9);
            %    com_vdw_new = com_vdw_new.^(1/0.9);
            %end
            C2arXjz_atom=[];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
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
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            
        elseif ligand(jz,4) == 6
            skhs1 = skhs1^(1/35);
            com_vdw_new = com_vdw_new.^(1/35);
            connect_num = connectdetect2(ligand(jz,:),O3_atom,R1);
            if connect_num>=1
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R1);
            if connect_num>=1
                
                skhs1 = skhs1^(1/1);
                com_vdw_new = com_vdw_new.^(1/1);
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,O2_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/1.4);
                com_vdw_new = com_vdw_new.^(1/1.4);
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/1.4);
                com_vdw_new = com_vdw_new.^(1/1.4);
            end
            
            C2arXjz_atom=[];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
        elseif ligand(jz,4) >= 7 && ligand(jz,4) <= 9
            skhs1 = skhs1^(1/28);
            com_vdw_new = com_vdw_new.^(1/28);
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/1.4);
                com_vdw_new = com_vdw_new.^(1/1.4);
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Nall_atom);
            if angle_mark>=1
                skhs1 = skhs1^(1/1.4);
                com_vdw_new = com_vdw_new.^(1/1.4);
            end
            
            angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,C2ar_atom);
            if angle_mark>=1
                
                skhs1 = skhs1^(1/2.1);
                com_vdw_new = com_vdw_new.^(1/2.1);
            end
            
            if ligand(jz,4) == 9
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            connect_num1 = connectdetect2(ligand(jz,:),Nar_atom,R1);
            angle_mark1 = angle_find_W(ligand(jz,:),Nar_atom,O2_atom);
            if angle_mark1>=1 && connect_num1>=1
                
                skhs1 = skhs1^(1/2);
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
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
            
            
            
            
            
            H_donor_atom_input = [];
            H_donor_atom_input = H_donor_atom(dist(H_donor_atom,ligand(jz,:))<3.0,:);
            intra_H_bond_atom = [];
            intra_H_bond_atom = H_donor_atom_input(dist(H_donor_atom_input,ligand(jz,:))<3.0 & dist(H_donor_atom_input,ligand(jz,:))>2.55,:);
            
            if ~isempty(intra_H_bond_atom)
                skhs1 = skhs1^(1/2);
                com_vdw_new = com_vdw_new.^(1/2);
            end
            %%%torsion_mark = torsion_find(ligand(jz,:),C2ar_atom,C2ar_atom);
            %%%if torsion_mark>=1
            %%%    skhs1 = skhs1^(1/tr2);
            %%%    com_vdw_new = com_vdw_new.^(1/tr2);
            %%%end
            
        elseif ligand(jz,4) == 20 || ligand(jz,4) == 21 || ligand(jz,4) == 14
            
            
            skhs1 = skhs1^(1/36);
            com_vdw_new = com_vdw_new.^(1/36);
            
            connect_num = connectdetect2(ligand(jz,:),Other_atom,R2);
            connect_numS = connectdetect2(ligand(jz,:),SP_atom,R2);
            if connect_numS >= 1
                skhs1 = skhs1^(1/1.4);
                com_vdw_new = com_vdw_new.^(1/1.4);
            elseif connect_numS == 0 && connect_num >= 2
                skhs1 = skhs1^(1/1);
                com_vdw_new = com_vdw_new.^(1/1);
            end
            
            
            
            
            
            
            
            
            
            
            %angle_mark = angle_find_W(ligand(jz,:),C2ar_atom,Oall_atom);
            
            
        elseif ligand(jz,4) == 17 || ligand(jz,4) == 18
            skhs1 = skhs1^(1/67);
            com_vdw_new = com_vdw_new.^(1/67);
            angle_mark = angle_find_W(ligand(jz,:),C3_atom,Halo_atom);
            if angle_mark>=2
                skhs1 = skhs1^(1/3.7);
                com_vdw_new = com_vdw_new.^(1/3.7);
            end
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R2);
            if connect_num>=1
                skhs1 = skhs1^(1/3.7);
                com_vdw_new = com_vdw_new.^(1/3.7);
            end
            
            C2arXjz_atom=[];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
        elseif ligand(jz,4) == 15
            skhs1 = skhs1^(-1/81);
            com_vdw_new = com_vdw_new.^(-1/81);
            angle_mark = angle_find_W(ligand(jz,:),C3_atom,Halo_atom);
            if angle_mark>=2
                skhs1 = skhs1^(1/3.7);
                com_vdw_new = com_vdw_new.^(1/3.7);
            end
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R2);
            if connect_num>=1
                skhs1 = skhs1^(1/3.7);
                com_vdw_new = com_vdw_new.^(1/3.7);
            end
            
            C2arXjz_atom=[];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
        elseif ligand(jz,4) == 16
            skhs1 = skhs1^(1/119);
            com_vdw_new = com_vdw_new.^(1/119);
            angle_mark = angle_find_W(ligand(jz,:),C3_atom,Halo_atom);
            if angle_mark>=2
                skhs1 = skhs1^(1/3.7);
                com_vdw_new = com_vdw_new.^(1/3.7);
            end
            connect_num = connectdetect2(ligand(jz,:),C2ar_atom,R2);
            if connect_num>=1
                skhs1 = skhs1^(1/3.7);
                com_vdw_new = com_vdw_new.^(1/3.7);
            end
            C2arXjz_atom=[];
            oppo_atom_final1 = [];
            oppo_atom_final2 = [];
            C2arXjz_atom = ligand(dist(ligand,ligand(jz,:))<=R1 & dist(ligand,ligand(jz,:))>0.5 & (ligand(:,4)==30 | ligand(:,4)==33),:);
            [oppo_atom_final1,oppo_atom_final2,torsion_mark] = oppo_atom_find(C2arXjz_atom,C2ar_atom,C2arX_atom,Polar_Aro_atom,Polar_Aro_atom);
            if ~isempty(oppo_atom_final1) && ~isempty(oppo_atom_final2) && oppo_atom_final1(1,4) == oppo_atom_final2(1,4) && oppo_atom_final1(1,7) == ligand(jz,7)
                skhs1 = skhs1^(1/1.5);
                com_vdw_new = com_vdw_new.^(1/1.5);
            end
            
        else
            skhs1 = skhs1^(1/-200);
            com_vdw_new = com_vdw_new.^(1/-200);
        end
        
        if ligand(jz,4) == 27
            skhs1 = skhs1^(1/-1.6);
            com_vdw_new = com_vdw_new.^(1/-1.6);
        end
        %if ligand(jz,4) == 30
        %    skhs1 = skhs1^(1/0.2);
        %    com_vdw_new = com_vdw_new.^(1/0.2);
        %end
        
        %%%cf_vdw_new=sum(cf_vdw);
        if ring_mark == 1 && ring_size == 3 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
            %if ring_mark == 1 && ring_size == 3 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
            skhs1 = skhs1^(1/2.8);
            com_vdw_new = com_vdw_new^(1/2.8);
            %skhs1 = skhs1^(1/3);
            %com_vdw_new = com_vdw_new^(1/3);
        elseif ring_mark == 1 && ring_size == 4 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
            %elseif ring_mark == 1 && ring_size == 4 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
            skhs1 = skhs1^(1/2.4);
            com_vdw_new = com_vdw_new^(1/2.4);
            %skhs1 = skhs1^(1/2.8);
            %com_vdw_new = com_vdw_new^(1/2.8);
        elseif ring_mark == 1 && ring_size == 5 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
            %elseif ring_mark == 1 && ring_size == 5 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
            skhs1 = skhs1^(1/2.2);
            com_vdw_new = com_vdw_new^(1/2.2);
            %skhs1 = skhs1^(1/2.6);
            %com_vdw_new = com_vdw_new^(1/2.6);
        elseif ring_mark == 1 && ring_size == 6 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || (ligand(jz,4)==33 & connectdetect2(ligand(jz,:),C2ar_atom,R1) > 0))
            %elseif ring_mark == 1 && ring_size == 6 && (ligand(jz,4)==1 || ligand(jz,4)==2 || ligand(jz,4)==4 || ligand(jz,4)==27 || ligand(jz,4)==30 || ligand(jz,4)==33)
            skhs1 = skhs1^(1/2);
            com_vdw_new = com_vdw_new^(1/2);
            %skhs1 = skhs1^(1/2);
            %com_vdw_new = com_vdw_new^(1/2);
        end
        ring_mark_ligand(jz,1)=ring_mark;
        ring_mark_ligand(jz,2)=ring_size;
        %if com_vdw_new>0
        Part_Matr_san=Part_Matr_san+log(skhs1)+log(com_vdw_new);
        %%%Freq_Matr_san=Freq_Matr_san+log(cf_vdw_new);
        %Exam1(jz,pe)=log(skhs1);
        %Exam2(jz,pe)=log(com_vdw_new);
        %pe=pe+1;
        skhs2 = skhs2^(1/19);
        com_vdw_new2 = com_vdw_new2.^(1/19);
        Part_Matr_san2=Part_Matr_san2+log(skhs2)+log(com_vdw_new2);
    end
    
    
end
SA_list_sum = sum(SA_list(:,1));
Part_Matr_san1 = -0.5918*Part_Matr_san;
Part_Matr_san2 = -0.5918*Part_Matr_san2;
%end
%end
%end
%Part_Matr_san
%%%vn(vn(:,1) == 0,1)=1;
%%%VN1 = sum(log(vn));

%max_san1=sum(sum(Freq_Matr.*Part_Matr_san));
%max_san1=Part_Matr_san;

%AdGs=-0.5918*(max_san1);

%out = VN1;
%if strcmp(ligand_pot,'0114NNd')==1
%    aaaaaaaaaaaaaaaa=1
%    pause
%end


%clear Part_Matr Freq_Matr Part_Matr_lig Freq_Matr_lig Part_Matr_com Freq_Matr_com Part_Matr_com Freq_Matr_com Part_Matr_ang Freq_Matr_ang Part_Matr_non Freq_Matr_non

%clear protein_name protein_atom_radius protein_atom_hb protein_acceptor_hb_angle protein_donor_hb_angle protein_atom_charge protein_backbone protein_chainID protein_resseq






