%%%%% Rading a pdb file %%%%%%%%%%%%%%%%%%%%%
tline=fgetl(fidin);
o=1;
while ~feof(fidin)
    if (strcmp(tline(1:3),'TER')==1)
        tline=fgetl(fidin);
    
    elseif ((strcmp(tline(1:4),'ATOM')==1)) 
        if (strfind(tline,'ILE')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
	elseif (strfind(tline,'HOH')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);

        elseif (strfind(tline,'AS')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
        elseif (strfind(tline,'SER')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);    

	elseif (strfind(tline,'ALA')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);

	elseif (strfind(tline,'PRO')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
	elseif (strfind(tline,'GL')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'VAL')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'THR')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'CY')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'HI')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin); 
        
    elseif (strfind(tline,'LEU')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'PHE')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'TRP')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'TYR')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'ARG')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'LYS')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
        
    elseif (strfind(tline,'MET')~=0)
        protein_atom(o,1)=str2num(tline(30:38));
        protein_atom(o,2)=str2num(tline(40:46));
        protein_atom(o,3)=str2num(tline(48:54));
        protein_res(o,:)=(tline(18:20));
        protein_resseq(o,:)=(tline(23:28));
        AtomName(o,:)=(tline(13:16));
        chainID(o,:) = (tline(22:22));
        o=o+1;
        tline=fgetl(fidin);
	end

    else
        tline=fgetl(fidin);
    end
end

  protein_atom1=protein_atom(:,1);
  protein_atom2=protein_atom(:,2);
  protein_atom3=protein_atom(:,3);
  fclose(fidin);
  clear protein_atom;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  clear pro atom seq chain h_bond res;
  k=1;
for i=1:length(protein_atom1)
    s=strtrim(AtomName(i,:));
    bb=strtrim(protein_res(i,:));
    if ((strncmp(s,'H',1)==1) || (strncmp(s,'1',1)==1) || (strncmp(s,'2',1)==1) || (strncmp(s,'3',1)==1)) %%&& (strncmp(bb,'HOH',3) ~= 1) 
        continue;
    else
        pro(k,:)= [protein_atom1(i,1) protein_atom2(i,1) protein_atom3(i,1)];
        atom11(k,:)=AtomName(i,:); res(k,:)=protein_res(i,:);
        seq(k,:)=protein_resseq(i,:);
        chain(k,:)=chainID(i,:);
        h_bond(k,1)=0;
        %%%%%% donors %%%
        if ((strncmp(s,'N',1)==1) && (length(s)==1)) || ((strncmp(s,'H1',2)==1) && (strncmp(bb,'HOH',3)==1)) || ((strncmp(s,'H2',2)==1) && (strncmp(bb,'HOH',3)==1)) || ((strncmp(s,'NZ',2)==1) && (strncmp(bb,'LYS',3)==1)) || ((strncmp(s,'ND1',3)==1) && (strncmp(bb,'HIS',3)==1)) || ((strncmp(s,'NE2',3)==1) && (strncmp(bb,'HIE',3)==1)) || ((strncmp(s,'NE1',3)==1) && (strncmp(bb,'TRP',3)==1)) || ((strncmp(s,'ND2',3)==1) && (strncmp(bb,'ASN',3)==1)) || ((strncmp(s,'NE2',3)==1) && (strncmp(bb,'GLN',3)==1)) || (((strncmp(s,'NH1',3)==1) || (strncmp(s,'NH2',3)==1) || (strncmp(s,'NE',2)==1)) && (strncmp(bb,'ARG',3)==1))
            h_bond(k,1)=1;
            %%%% acceptor %%%%%
        elseif ((strncmp(s,'O',1)==1) && (length(s)==1)) || ((strncmp(s,'O',1)==1) && (strncmp(bb,'HOH',3)==1)) || (((strncmp(s,'OD1',3)==1) || (strncmp(s,'OD2',3)==1)) && ((strncmp(bb,'ASN',3)==1) || (strncmp(bb,'ASP',3)==1))) || (((strncmp(s,'OE1',3)==1) || (strncmp(s,'OE2',3)==1)) && ((strncmp(bb,'GLN',3)==1) || (strncmp(bb,'GLU',3)==1))) || ((strncmp(s,'SD',2)==1) && (strncmp(bb,'MET',3)==1))  || ((strncmp(s,'ND1',3)==1) && (strncmp(bb,'HIE',3)==1) )
            h_bond(k,1)=2;
            %%%% both %%%%%%%
        elseif ((strncmp(s,'SG',2)==1) && (strncmp(bb,'CYS',3)==1)) || ((strncmp(s,'OH',2)==1) &&  (strncmp(bb,'TYR',3)==1)) || ((strncmp(s,'OG',2)==1) && (strncmp(bb,'SER',3)==1)) || ((strncmp(s,'OG1',3)==1) && (strncmp(bb,'THR',3)==1))  %|| ((strncmp(s,'NE2',3)==1) && (strncmp(bb,'HIS',3)==1) )
            h_bond(k,1)=3;
        else h_bond(k,1)=0;
        end
        
        k=k+1;
    end
end
