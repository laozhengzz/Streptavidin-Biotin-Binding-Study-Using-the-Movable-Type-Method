%%%%% Read a ligand mol2 file  %%%%%%%%%%%%%%%

%clear atom_conj_num atom name name_str ligand_atom_num atom_contact contact_type contact_name atom_contact_num mass;
fidin=fopen(ligand_dir(1,:));
tline=fgetl(fidin);
k=1;

while ~feof(fidin)
    if length(tline)>=70
        latom(k,1)=str2num(tline(18:26));
        latom(k,2)=str2num(tline(28:36));
        latom(k,3)=str2num(tline(38:46));
        name_str(k,1:5)=tline(48:52);
        %ligand_atom_num(k,1)=str2num(tline(1:7));
        k=k+1;
        tline=fgetl(fidin);
%     elseif (length(tline)==13)&(strcmp(tline(1:13), '@<TRIPOS>BOND')==1)
%         tline=fgetl(fidin);
%         while (strcmp(tline, '@<TRIPOS>SUBSTRUCTURE')~=1)
% %            tline=fgetl(fidin);
%             atom_contact(l,1)=str2num(tline(9:11));
%             l=l+1;
%             atom_contact(l,1)=str2num(tline(14:16));
%             l=l+1;
%             contact_type(m,1)=str2num(tline(9:11));
%             contact_type(m,2)=str2num(tline(14:16));
%             contact_name(m,1:2)=tline(18:19);
%             m=m+1;
%             tline=fgetl(fidin);
%             if feof(fidin)==1
%                 break
%             end
%         end
    elseif (length(tline)==21)&(strcmp(tline(1:21), '@<TRIPOS>SUBSTRUCTURE')==1)
        break
    else
        tline=fgetl(fidin);
    end
end

size_atom=size(latom);


fclose(fidin);

%ligand initialization
        

clear name;  

size_name_str=size(name_str);
name=zeros(length(name_str),1);

for i=1:1:size_name_str(1,1)
    if strcmp(name_str(i,1:5), 'H    ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'C.3  ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'C.2  ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'C.1  ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'C.ar ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'O.3  ')==1
        name(i,1)=3;
    elseif strcmp(name_str(i,1:5), 'O.2  ')==1
        name(i,1)=2;   
    elseif strcmp(name_str(i,1:5), 'N.3  ')==1
        name(i,1)=3;
    elseif strcmp(name_str(i,1:5), 'N.2  ')==1
        name(i,1)=1;
    elseif strcmp(name_str(i,1:5), 'N.1  ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'N.ar ')==1
        name(i,1)=1;
    elseif strcmp(name_str(i,1:5), 'N.am ')==1
        name(i,1)=1;
    elseif strcmp(name_str(i,1:5), 'N.pl3')==1
        name(i,1)=1;
    elseif strcmp(name_str(i,1:5), 'S.3  ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'S.2  ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'S.o  ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'S.o2 ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'F    ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'Cl   ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'Br   ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'I    ')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'P.3  ')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'C.cat')==1
        name(i,1)=0;
    elseif strcmp(name_str(i,1:5), 'O.co2')==1
        name(i,1)=2;
    elseif strcmp(name_str(i,1:5), 'N.4  ')==1
        name(i,1)=1;
    end
end


 %clear pro atom seq chain h_bond res;
  k=1;
for i=1:length(latom)
    s=strtrim(name_str(i,:));
    %bb=strtrim(protein_res(i,:));
    if ((strncmp(s,'H',1)==1) || (strncmp(s,'1',1)==1) || (strncmp(s,'2',1)==1) || (strncmp(s,'3',1)==1)) 
        continue;
    else
        lig_atom(k,:)= [latom(i,1) latom(i,2) latom(i,3)];
        lig_hbond(k,1)=name(i,1);
        k=k+1;
    end
end
