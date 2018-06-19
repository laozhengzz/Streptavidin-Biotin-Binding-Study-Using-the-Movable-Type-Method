%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Mtflex for backbone and side chain flexibility                    
%  Written by: Nupur Bansal, Merz Group                  
%  Dated: 10/06/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear; clc;
function mtflex_parallel_read(protein_dir,pocket_dir,ligand_dir)

% %Begin=29501; Constant=500;
% load tor.mat;
% 
% PR=P_all(:,Begin:Begin+Constant-1);
% EN=E_all(:,Begin:Begin+Constant-1);
% TO=Tor_all(:,Begin:Begin+Constant-1);

cut=0.11;
%delete pos_*.pdb *.txt;
%% Read in your protein pdb and also the pocket pdb

%protein_dir=('/Users/Nupur/Downloads/1iep_protein.pdb');
%pocket_dir=('/Users/Nupur/Downloads/1iep_pocket.pdb');
%ligand_dir=('/Users/Nupur/Downloads/1iep_ligand.mol2');
% %protein_dir=('/Users/Nupur/Receptor_flex/4OBO-monomer_Stage1.pdb');
% %pocket_dir=('/Users/Nupur/Receptor_flex/4obo_pocket.pdb');
%  protein_dir=('/Users/Nupur/apo_homo_conformer/1mk5_1swc_try/1swc_protein.pdb');
%  pocket_dir=('/Users/Nupur/apo_homo_conformer/1mk5_1swc_try/1swc_pocket.pdb');

%% Open your protein & pocket pdb and read out the information from the file
%%%%  *Warning: Only heavy atoms of 20-Amino Acids will be read*


%% Protein
fidin = fopen (protein_dir);
pdb_read;
pro_coor=pro; pro_res=res; pro_seq=seq; pro_atom=atom11; pro_hbond=h_bond;
pro_chain=chain;
size_protein=size(pro_coor);


%% Pocket
fidin = fopen (pocket_dir);
pdb_read;
poc_coor=pro; poc_res=res; poc_seq=seq; poc_atom=atom11; poc_hbond=h_bond;
poc_chain=chain;
size_pocket=size(poc_coor);

%% ligand
lig_read;

% %%
% %%%%%%%%% scan disulphide bond in crystal str. %%%%%%%%%%%%%%%%%%%%
% disulphide_bond=zeros(size_protein(1,1),1);
% 
% for i = 1:size_protein(1,1)
%     s=strtrim(pro_atom(i,:));
%     bb=strtrim(pro_res(i,:));
%     for j = 1:size_protein(1,1)
%         s1=strtrim(pro_atom(j,:));
%         bb1=strtrim(pro_res(j,:));
%         if j ~= i
%             if (((strncmp(s,'SG',2)==1) && (strncmp(bb,'CYS',3)==1)) ...
%                     && ((strncmp(s1,'SG',2)==1) && (strncmp(bb1,'CYS',3)==1))) ...
%             && (sqrt((pro_coor(j,1)-pro_coor(i,1))^2 +(pro_coor(j,2)-pro_coor(i,2))^2+(pro_coor(j,3)-pro_coor(i,3))^2)<=2.2)
%                 disulphide_bond(i)=disulphide_bond(i)+1;
%                 break;
%             end
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%% storing the rest of atoms as rest matrix %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% declaring empty arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr=zeros(size_pocket(1,1),1);
name_res=char(zeros(size_pocket(1,1),3));
seq_res=char(zeros(size_pocket(1,1),6));
chain_res=char(zeros(size_pocket(1,1),1));
adj_p=char(zeros(size_pocket(1,1),6));
adj_n=adj_p; adj_pc=char(zeros(size_pocket(1,1),1));
adj_nc=adj_pc;
[k1,k2]=size(pro_coor);
num_always=k2/3;
back_n=zeros(size_pocket(1,1),k2);
back_ca=zeros(size_pocket(1,1),k2);
back_c=zeros(size_pocket(1,1),k2);
back_o=zeros(size_pocket(1,1),k2);
back_cb=zeros(size_pocket(1,1),k2);

rest_coor=zeros((size_protein(1,1)-size_pocket(1,1)),3);
rest_res=char(zeros((size_protein(1,1)-size_pocket(1,1)),3));
rest_atom=char(zeros((size_protein(1,1)-size_pocket(1,1)),4));
rest_chain=char(zeros((size_protein(1,1)-size_pocket(1,1)),1));
rest_seq=char(zeros((size_protein(1,1)-size_pocket(1,1)),6));
rest_hbond=zeros((size_protein(1,1)-size_pocket(1,1)),1);
rest_num=zeros((size_protein(1,1)-size_pocket(1,1)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear c1 c2;
clear i j;
a=0; b=0; c=0;
k=1;

for j = 1:size_protein(1,1)
    ff=0;
    
    for i = 1:size_pocket(1,1)
        s1=strtrim(poc_atom(i,:));
        bb1=strtrim(poc_res(i,:));
        if (poc_coor(i,1)==pro_coor(j,1) && (poc_coor(i,2)==pro_coor(j,2)) && (poc_coor(i,3)==pro_coor(j,3)))
            ff=1;
            %
            
	%%% to note the name of the residue & keeping backbone atoms %%%
            if (strncmp(s1,'N',1)==1) && (length(s1)==1)
                a=a+1;

                adj_n(a,:)=pro_seq(j-1,:); adj_nc(a,:)=pro_chain(j-1,:);
                chain_res(a,:)=pro_chain(j-1,:);
                

                [n_atoms, name_atoms, d1_atoms, d2_atoms, d3_atoms, hb] = aa_info(bb1);
                nr(a,1)=n_atoms;


                adj_p(a,:)=pro_seq(j+n_atoms,:);
                adj_pc(a,:)=pro_chain(j+n_atoms,:);

                name_res(a,:) = poc_res(i,:);
                seq_res(a,:)  = poc_seq(i,:);
            end

        end
    end

    if ff==0
        b=b+1;
        rest_coor(b,:)=pro_coor(j,:); rest_res(b,:)=pro_res(j,:); rest_atom(b,:)=pro_atom(j,:);
        rest_chain(b,:)=pro_chain(j,:);rest_seq(b,:)=pro_seq(j,:); rest_hbond(b,:)=pro_hbond(j,:);
        rest_num(b,1)=j;
    end
end
name_res = name_res(1:a,:);
seq_res = seq_res(1:a,:);
chain_res = chain_res(1:a,:);
adj_n = adj_n(1:a,:); adj_p = adj_p(1:a,:);
adj_nc = adj_nc(1:a,:); adj_pc = adj_pc(1:a,:);
num_res=a;
nr = nr(1:a,1);
disp(num_res);

%%%% Declare variables to fill in backbone information
%%%% for adding coordinates
%%%   We will fill the backbone coordinates of N,Ca,C and O in this matrix.
Point_back_1 = zeros(4,3);
Point_back_n = zeros(4,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% segregate backbone coordinates of the i-1 residue 
%%%%% where i is the first residue of the loop
i=1;   
a=strtrim(adj_n(i,:));
c=strtrim(adj_nc(i,:));

for j=1:size_protein(1,1)
    s1=strtrim(pro_atom(j,:));
    a1=strtrim(pro_seq(j,:));
    c1=strtrim(pro_chain(j,:));
    
    if ((strncmp(s1,'N',1)==1) && (length(s1)==1)) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_1(1,1:3)=pro_coor(j,:);
        
    elseif ((strncmp(s1,'C',1)==1) && (length(s1)==1)) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_1(3,1:3)=pro_coor(j,:);
        
    elseif ((strncmp(s1,'O',1)==1) && (length(s1)==1)) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_1(4,1:3)=pro_coor(j,:);
        
    elseif (strncmp(s1,'CA',2)==1) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_1(2,1:3)=pro_coor(j,:);         
    end    
end

%%%%% segregate backbone coordinates of the n+1 residue 
%%%%% where n is the last residue of the loop

i=num_res;   
a=strtrim(adj_p(i,:));
c=strtrim(adj_pc(i,:));

for j=1:size_protein(1,1)
    s1=strtrim(pro_atom(j,:));
    a1=strtrim(pro_seq(j,:));
    c1=strtrim(pro_chain(j,:));
    
    if ((strncmp(s1,'N',1)==1) && (length(s1)==1)) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_n(1,1:3)=pro_coor(j,:);
        
    elseif ((strncmp(s1,'C',1)==1) && (length(s1)==1)) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_n(3,1:3)=pro_coor(j,:);
        
    elseif ((strncmp(s1,'O',1)==1) && (length(s1)==1)) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_n(4,1:3)=pro_coor(j,:);
        
    elseif (strncmp(s1,'CA',2)==1) && (strncmp(a,a1,3)==1) && (strncmp(c,c1,3)==1)
        Point_back_n(2,1:3)=pro_coor(j,:);         
    end    
end




clear a bb ff;

save chk.mat;
end
