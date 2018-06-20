function out = GARF_MT_EXE_copy(home_dir)
%clear all
warning off


load torsion_AF1
load torsion_list_final_AF
torsion_list_num_AF = contact_list_final_AF_format(torsion_list_final_AF);
load vdw_d_all_AF
load vdw_all_AF_contact_list_num

load nonb_d_f_refined
GARF_Potential
%run GARF_nonbonding_PP

Solvation_Mat
%run GARF_inter_PP

%run LL_torsion_FF

rd = 0.005;

R_srch=0.2;
R_L_tor_srch=0.1;
hnstate = 0.2/rd;
halva=hnstate;
halv=100;
halv1=10;

const=halva*2;
const1=halv1*2;
column_num=100;
size_nonb_d_f_refined=size(nonb_d_f_refined);
for i=2:size_nonb_d_f_refined(1,2)
    nonb_d_f_refined(nonb_d_f_refined(:,i)<=10^0,i)=10^0;
end
size_vdw_d=size(vdw_d);
for i=2:size_vdw_d(1,2)
    vdw_d(vdw_d(:,i)<=10^-10,i)=10^-10;
end

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



carbon='C';
oxygen='O';
nitrogen='N';
sulfer='S';
fluorin='F';
chlorine='Cl';
bromine='Br';
iodine='I';
phosphor='P';
hydrogen='H';

RcutoffPP = 6;
RcutoffPL = 6;

count=1

%APO_protein_dir1 = 'C:\Users\John Zheng\Dropbox\Test_cases\Streptavidin_apo_holo\Apo_protein\';
%list_ApoP=dir(HOLO_protein_dir1);

%home_dir = 'F:\Streptavidin_apo_holo\crystal_close\';
list=dir(home_dir);
size_list=size(list);
%Aout=zeros(200,6);
%Bout=zeros(200,28);

tn=0;
clear test_name
%for zh=15:15
fp1=fopen(strcat(home_dir,'PL_MT_tet_binding_0.15.txt'),'wt');
fp2=fopen(strcat(home_dir,'P_MT_tet_0.15.txt'),'wt');
fp3=fopen(strcat(home_dir,'PL_MT_tet_0.15.txt'),'wt');
%fp4=fopen(strcat(home_dir,'solv.txt'),'wt'); 

for zh=3:size_list(1,1)
    
    tic
    if length(list(zh,1).name)<8
        continue
    elseif length(list(zh,1).name)>=8
        list_name=list(zh,1).name;
        %list_name{zh,1}=list_name;
        %if strcmp(list_name(end-7:end),'holo.pdb')==1
        if strcmp(list_name(end-3:end),'.pdb')==1 
	   protein_dir=[ strcat(home_dir,list_name)]
	   [protein,protein_Namelist,protein_sol] = protein_define_final_sol(protein_dir);

            %[Apo_protein,Apo_protein_Namelist,Apo_protein_sol] = protein_define_final_sol(Apo_protein_dir);
           Apo_protein=protein; Apo_protein_Namelist=protein_Namelist;
           Apo_protein_sol=protein_sol;
            protein(protein(:,5)==121,:)=[];
            %Apo_protein(Apo_protein(:,5)==121,:)=[];
            Apo_protein=protein;
          
	  for zh1=1:5
	    t1=sprintf('_btn_%d.mol2',zh1);  
 
            id_name = list_name(1:end-4);
            ligand_pot = strcat(id_name,t1);
%%%%%%%	   ligand_pot='biotin.mol2'; 
           ligand_dir=[ strcat(home_dir,ligand_pot)]
            
            %Apo_protein_pot = strcat(id_name,'_apo.pdb');
            %Apo_protein_dir=[ strcat(home_dir,Apo_protein_pot)];
            Apo_protein_dir=protein_dir;
            [ligand,name_str] = ligand_define_sol(ligand_dir);
            
            ring_atom = ring_judge(ligand);
            
            complex = [protein_sol
                ligand];
            
            [ZE_P_san,ZE_P_water,SA_P_sum] = KMTISM_function_new(protein_sol);
            toc
            tic
            %[ZE_AP_san,ZE_AP_water,SA_AP_sum] = KMTISM_function_new(Apo_protein_sol);
            ZE_AP_san=ZE_P_san; ZE_AP_water=ZE_P_water; SA_AP_sum=SA_P_sum;
	     toc
            tic
            [ZE_L_san,ZE_L_water,SA_L_sum] = KMTISM_function_new(ligand);
            toc
            tic
            [ZE_PL_san,ZE_PL_water,SA_PL_sum] = KMTISM_function_new(complex);
            toc
            
            
            tic
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%protein%%%%%%%%%%%%%%%%%%%%%%%%%%
            angle_list_P = angle_find(protein);
            [torsion_list_P_Namelist,torsion_list_P_all] = torsion_find_protein_final_AFA(protein,protein_Namelist);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Apo_protein%%%%%%%%%%%%%%%%%%%%%%
           % angle_list_AP = angle_find(Apo_protein);
           angle_list_AP=angle_list_P;
	  % [torsion_list_AP_Namelist,torsion_list_AP_all] = torsion_find_protein_final_AFA(Apo_protein,Apo_protein_Namelist);
            torsion_list_AP_Namelist=torsion_list_P_Namelist;
	    torsion_list_AP_all=torsion_list_P_all;
	    [Part_Matr_PTor,CNNPTor,VP,EPT_PTor,Part_Matr_APTor,CNNAPTor,VAP,EPT_APTor] = torsion_potential_AFA(torsion_list_P_all,torsion_list_P_Namelist,torsion_list_AP_all,torsion_list_AP_Namelist,torsion_list_num_AF,torsion_AF1,R_L_tor_srch);
            
            toc
            tic
            
            
            [Part_Matr_PNonB1,CNNPNonB1,VPNonB1]  = intra_nonbond_potential_pp_GARFv1(protein,protein_Namelist,torsion_list_P_all,angle_list_P,vdw_d_all_AF,vdw_all_AF_contact_list_num,nonb_d_f_refined,RcutoffPP);
            %Part_Matr_PNonB2 = Part_Matr_PNonB2(1,1);
            toc
            tic
            
            %[Part_Matr_APNonB1,CNNAPNonB1,VAPNonB1]  = intra_nonbond_potential_pp_GARFv1(Apo_protein,Apo_protein_Namelist,torsion_list_AP_all,angle_list_AP,vdw_d_all_AF,vdw_all_AF_contact_list_num,nonb_d_f_refined,RcutoffPP);
            %Part_Matr_APNonB2 = Part_Matr_APNonB2(1,1);
	    Part_Matr_APNonB1=Part_Matr_PNonB1; CNNAPNonB1=CNNPNonB1;
	    VAPNonB1=VPNonB1;
            toc
            
            tic
            [protein_starterA,protein_starterA_Namelist,protein_starterB] = pocket2find_PL_AFA(protein,protein_Namelist,ligand,RcutoffPL);
            [Part_Matr_com1,Part_Matr_com2,CNN1,CNN2,VN1,VN2] = inter_potential_PL_GARF(protein_starterA,protein_starterB,vdw_d);
            
            ZE_water = ZE_PL_water - ZE_P_water - ZE_L_water;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            ZE_part_com1 = -0.5918*(Part_Matr_com1 - VN1);
            ZN_part_com1 = ((sqrt(CNN1*8+1)+1)/2-4)*log(2)+((sqrt(CNN1*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1*8+1)+1)/2-1)*log(const);
            
            ZE_part_com2 = -0.5918*(Part_Matr_com2 - VN2);
            ZN_part_com2 = ((sqrt(CNN2*8+1)+1)/2-4)*log(2)+((sqrt(CNN2*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN2*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN2*8+1)+1)/2-1)*log(const);
            
            ZE_part_PTor = -0.5918*(Part_Matr_PTor - VP);
            ZN_part_PTor = ((sqrt(CNNPTor*8+1)+1)/2-4)*log(2)+((sqrt(CNNPTor*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNNPTor*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNNPTor*8+1)+1)/2-1)*log(const);
            
            ZE_part_PNonB1 = -0.5918*(Part_Matr_PNonB1 - VPNonB1);
            ZN_part_PNonB1 = ((sqrt(CNNPNonB1*8+1)+1)/2-4)*log(2)+((sqrt(CNNPNonB1*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNNPNonB1*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNNPNonB1*8+1)+1)/2-1)*log(const);
            
            ZE_part_ApoPTor = -0.5918*(Part_Matr_APTor - VAP);
            ZN_part_ApoPTor = ((sqrt(CNNAPTor*8+1)+1)/2-4)*log(2)+((sqrt(CNNAPTor*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNNAPTor*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNNAPTor*8+1)+1)/2-1)*log(const);
            
            ZE_part_ApoPNonB1 = -0.5918*(Part_Matr_APNonB1 - VAPNonB1);
            ZN_part_ApoPNonB1 = ((sqrt(CNNAPNonB1*8+1)+1)/2-4)*log(2)+((sqrt(CNNAPNonB1*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNNAPNonB1*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNNAPNonB1*8+1)+1)/2-1)*log(const);
            
            ZE_sol = ZE_PL_san - ZE_AP_san - ZE_L_san;
            %fprintf(fp4,' % s  % 5d\n',protein_dir,ZE_sol);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Protein_energy(count,1) =  ZE_part_PTor + ZN_part_PTor*-0.5918 + (ZE_part_PNonB1 + ZN_part_PNonB1*-0.5918)/(RcutoffPP*10) + ZE_P_san;
        	Ptor=ZE_part_PTor + (ZN_part_PTor*-0.5918) ;
		PNon=ZE_part_PNonB1 + (ZN_part_PNonB1*-0.5918);  
	 % Protein_energy(count,2) =  ZE_part_ApoPTor + ZN_part_ApoPTor*-0.5918 + (ZE_part_ApoPNonB1 + ZN_part_ApoPNonB1*-0.5918)/(RcutoffPP*10) + ZE_AP_san-0.43*ZE_P_water;
           % Protein_energy(count,3) =  Protein_energy(count,1) - Protein_energy(count,2);
	 fprintf(fp2,' % s  % 5d  % 5d  % 5d  % 5d  % 5d  % 5d  % 5d\n',ligand_dir,Ptor,PNon,ZE_P_san,ZE_L_san,ZE_P_water,ZE_L_water,Protein_energy(count,1));
	
   ProteinLigand_energy(count,1) = (ZE_part_com1 + ZN_part_com1*-0.5918)/(RcutoffPL) + ZE_part_ApoPTor + ZN_part_ApoPTor*-0.5918 + (ZE_part_ApoPNonB1 + ZN_part_ApoPNonB1*-0.5918)/(RcutoffPP*10) + ZE_PL_san;
           com1=ZE_part_com1 + (ZN_part_com1*-0.5918); 
            
           fprintf(fp3,' % s  % 5d  % 5d  % 5d  % 5d  % 5d  % 5d\n',ligand_dir,Ptor,PNon,com1,ZE_PL_san,ZE_PL_water,ProteinLigand_energy(count,1)); 
            Aout(count,1)= (ZE_part_com1 + ZN_part_com1*-0.5918)/(RcutoffPL*10) + (ZE_part_PTor + ZN_part_PTor*-0.5918) - (ZE_part_ApoPTor + ZN_part_ApoPTor*-0.5918) + ((ZE_part_PNonB1 + ZN_part_PNonB1*-0.5918) - (ZE_part_ApoPNonB1 + ZN_part_ApoPNonB1*-0.5918))/(RcutoffPP*10) + ZE_sol -0.43*ZE_water;
            Aout(count,2)= (ZE_part_com1 + ZN_part_com1*-0.5918)/(RcutoffPL*10) + ZE_sol -0.1*ZE_water;
	   fprintf(fp1,' % s  % 5d  % 5d  % 5d  % 5d  % 5d\n',ligand_dir,com1,ZE_sol,ZE_water,Aout(count,1),Aout(count,2));            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %fprintf(fp,' % s  % 5d  % 5d\n',ligand_pot,Aout(count,1), Bout(count,1));
            count=count+1
            
	end            
        end
    end
end
%fclose(fp)

fclose(fp1);
fclose(fp2);
fclose(fp3);
%fclose(fp4);

