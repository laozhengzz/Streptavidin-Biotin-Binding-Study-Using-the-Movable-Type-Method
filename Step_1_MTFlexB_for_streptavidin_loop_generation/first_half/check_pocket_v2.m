function  out1 = check_pocket_v2(h_bond,atom_new,res_new,seq_new,chain,i,drume,l,y1)



i_now=i; drume_now=drume; l_now=l; y1_now=y1;atom_new_now=atom_new; res_new_now=res_new; 
seq_new_now=seq_new; 
chain_now=chain; h_bond_now=h_bond;

load chk.mat;
% t=sprintf('%d.mat',i);
% load t;

i=i_now; drume=drume_now; l=l_now; y1=y1_now; atom_new=atom_new_now; res_new=res_new_now; 
seq_new=seq_new_now; chain=chain_now; h_bond=h_bond_now;

clear i_now drume_now k_now l_now y1_now atom_new_now res_new_now seq_now chain_now h_bond_now;
dist_h=2.2; dist_nb=2.80;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Vander waal collision within pocket %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for f=1:j2
%j1=3*(f-1);
dummy=drume;
a=0; k1=1;
for x=y1+1:l
    s=strtrim(atom_new(x,:));
    bb=strtrim(res_new(x,:));
    u=strtrim(seq_new(x,:));
    c=strtrim(chain(x,:));
    k1=i;  numn=nr(i,1); u_1=strtrim(adj_p(i,:)); u_11=strtrim(adj_n(i,:)); c_1=strtrim(adj_pc(i,:)); c_11=strtrim(adj_nc(i,:));
    if ((strncmp(s,'H1',2)==1) || (strncmp(s,'H2',2)==1)) && (strncmp(bb,'HOH',3)==1)
        dist_h = 1.4;
        dist_nb = 2.37;
    end
    d1=dummy(x,1:3);
    mat_dist=pdist2(dummy,d1);
    %aa1p=pro_new(x,1:3);
    if mean(d1)==0
        continue;
    end
    h1=h_bond(x,1);
    for g=1:l
        if mat_dist(g,1) < 2.8
            if x~=g
                s1=strtrim(atom_new(g,:));
                bb1=strtrim(res_new(g,:));
                u1=strtrim(seq_new(g,:));
                c1=strtrim(chain(g,:));
                %aa2=pro_new(g,1:3);
                d2=dummy(g,1:3);
                
                if (strncmp(u_1,u1,4)==1) & (strncmp(c_1,c1,3)==1)
                    if (strncmp(bb1,'PRO',3)==1)
                        if ((strncmp(s,'N',1)==1) && (length(s)==1)) && ((strncmp(s1,'N',1)==1) && (length(s1)==1))
                            continue;
                        elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && (((strncmp(s1,'N',1)==1) && (length(s1)==1)) || (strncmp(s1,'CG',2)==1) || (strncmp(s1,'CA',2)==1) || (strncmp(s1,'CB',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)))
                            continue;
                        elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && (strncmp(s1,'CD',2)==1)
                            aa3=sqrt((d1(1,1)-d2(1,1))^2 + (d1(1,2)-d2(1,2))^2 + (d1(1,3)-d2(1,3))^2);
                            if aa3 >= 2.00
                                continue;
                            else a=1; break;
                            end
                        elseif ((strncmp(s,'O',1)==1) && (length(s)==1)) && ((strncmp(s1,'CA',2)==1) || ((strncmp(s1,'N',1)==1) && (length(s1)==1))|| (strncmp(s1,'CD',2)==1))
                            continue;
                        elseif (strncmp(s,'CA',2)==1) && ((strncmp(s1,'CD',2)==1) || (strncmp(s1,'CA',2)==1) || ((strncmp(s1,'N',1)==1) && (length(s1)==1)))
                            continue;
                        elseif (strncmp(s,'CB',2)==1) && ((strncmp(s1,'N',1)==1) && (length(s1)==1))
                            continue;
                        end
                        
                    else
                        if ((strncmp(s,'N',1)==1) && (length(s)==1)) && ((strncmp(s1,'N',1)==1) && (length(s1)==1))
                            continue;
                        elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && (((strncmp(s1,'N',1)==1) && (length(s1)==1)) || (strncmp(s1,'CA',2)==1) || (strncmp(s1,'CB',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)))
                            continue;
                        elseif ((strncmp(s,'O',1)==1) && (length(s)==1)) && ((strncmp(s1,'CA',2)==1) || ((strncmp(s1,'N',1)==1) && (length(s1)==1)))
                            continue;
                        elseif (strncmp(s,'CA',2)==1) && ((strncmp(s1,'CA',2)==1) || ((strncmp(s1,'N',1)==1) && (length(s1)==1)))
                            continue;
                        elseif (strncmp(s,'CB',2)==1) && ((strncmp(s1,'N',1)==1) && (length(s1)==1))
                            continue;
                        end
                    end
                    
                elseif (strncmp(u_11,u1,4)==1) & (strncmp(c_11,c1,3)==1)
                    if (strncmp(bb1,'PRO',3)==1)
                        if (((strncmp(s,'N',1)==1) && (length(s)==1))) && (((strncmp(s1,'N',1)==1) && (length(s1)==1)) || (strncmp(s1,'CA',2)==1) || (strncmp(s1,'CB',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)) || ((strncmp(s1,'O',1)==1) && (length(s1)==1)))
                            continue;
                        elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && ((strncmp(s1,'C',1)==1) && (length(s1)==1))
                            continue;
                        elseif (strncmp(s,'CA',2)==1) && ((strncmp(s1,'CA',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)) || ((strncmp(s1,'O',1)==1) && (length(s1)==1)))
                            continue;
                        elseif (strncmp(s,'CD',2)==1) && ((strncmp(s1,'CA',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)) || ((strncmp(s1,'O',1)==1) && (length(s1)==1)))
                            continue;
                        elseif (strncmp(s,'CB',2)==1) && ((strncmp(s1,'C',1)==1) && (length(s1)==1))
                            continue;
                        elseif (strncmp(s,'CG',2)==1) && ((strncmp(s1,'C',1)==1) && (length(s1)==1))
                            continue;
                        end
                        
                    else
                        if (((strncmp(s,'N',1)==1) && (length(s)==1))) && (((strncmp(s1,'N',1)==1) && (length(s1)==1)) || (strncmp(s1,'CA',2)==1) || (strncmp(s1,'CB',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)) || ((strncmp(s1,'O',1)==1) && (length(s1)==1)))
                            continue;
                        elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && ((strncmp(s1,'C',1)==1) && (length(s1)==1))
                            continue;
                        elseif (strncmp(s,'CA',2)==1) && ((strncmp(s1,'CA',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)) || ((strncmp(s1,'O',1)==1) && (length(s1)==1)))
                            continue;
                        elseif (strncmp(s,'CB',2)==1) && ((strncmp(s1,'C',1)==1) && (length(s1)==1))
                            continue;
                        end
                    end
                end
                % put proper restraints on N, Ca, Cb
                
                if ((strncmp(s,'N',1)==1) & (length(s)==1)) & (((strncmp(u,u1,4)==1) & (bb1==bb) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | ((strncmp(s1,'O',1)==1) & (length(s1)==1)) | (strncmp(s1,'CG',2)==1) | (strncmp(s1,'CG1',3)==1) | (strncmp(s1,'CG2',3)==1) | ((strncmp(s1,'CD',2)==1) & (strncmp(bb,'PRO',3)==1)) | ((strncmp(s1,'SG',2)==1) & (strncmp(bb,'CYS',3)==1)) | ((strncmp(s1,'OG',2)==1) & (strncmp(bb,'SER',3)==1)) | ((strncmp(s1,'OG1',3)==1) & (strncmp(bb,'THR',3)==1)))))
                    continue;
                elseif (strncmp(s,'CA',2)==1) & (((strncmp(u,u1,4)==1) & (bb1==bb) & (((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'N',1)==1) & (length(s1)==1)) | ((strncmp(s1,'O',1)==1) & (length(s1)==1)) | (strncmp(s1,'CG',2)==1) | (strncmp(s1,'CG1',3)==1) | (strncmp(s1,'CG2',3)==1) |(strncmp(s1,'SD',2)==1) | (strncmp(s1,'CD',2)==1) | (strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (((strncmp(s1,'OD1',3)==1) | (strncmp(s1,'OD2',3)==1) | (strncmp(s1,'ND2',3)==1)) & ((strncmp(bb,'ASP',3)==1) | (strncmp(bb,'ASN',3)==1))) | ((strncmp(s1,'OG',2)==1) & (strncmp(bb,'SER',3)==1)) | (((strncmp(s1,'ND1',3)==1) | (strncmp(s1,'CD2',3)==1)) & (strncmp(bb,'HIS',3)==1)) | ((strncmp(s1,'OG1',3)==1) & (strncmp(bb,'THR',3)==1)) | ((strncmp(s1,'SG',2)==1) & (strncmp(bb,'CYS',3)==1))  )))
                    continue;
                elseif (strncmp(s,'CB',2)==1) & (((strncmp(u,u1,4)==1) & (bb1==bb) & (((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | ((strncmp(s1,'N',1)==1) & (length(s1)==1)) | ((strncmp(s1,'O',1)==1) & (length(s1)==1)) | (strncmp(s1,'CG',2)==1) | (strncmp(s1,'CG1',3)==1) | (strncmp(s1,'CG2',3)==1) |(strncmp(s1,'SD',2)==1) | (strncmp(s1,'CD',2)==1) | (strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE',2)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE3',3)==1) | (((strncmp(s1,'OD1',3)==1) | (strncmp(s1,'OD2',3)==1) | (strncmp(s1,'ND2',3)==1) | (strncmp(s1,'NE2',3)==1) | (strncmp(s1,'OE1',3)==1) | (strncmp(s1,'OE2',3)==1)) & ((strncmp(bb,'ASP',3)==1) | (strncmp(bb,'ASN',3)==1) | (strncmp(bb,'GLN',3)==1) | (strncmp(bb,'GLU',3)==1))) | ((strncmp(s1,'OG',2)==1) & (strncmp(bb,'SER',3)==1)) | (((strncmp(s1,'ND1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'NE2',3)==1)) & (strncmp(bb,'HIS',3)==1)) | ((strncmp(s1,'OG1',3)==1) & (strncmp(bb,'THR',3)==1)) | ((strncmp(s1,'NE',2)==1) & (strncmp(bb,'ARG',3)==1)) | ((strncmp(s1,'SG',2)==1) & (strncmp(bb,'CYS',3)==1)) | ((strncmp(s1,'NE1',3)==1) & (strncmp(bb,'TRP',3)==1)) )))
                    continue;
                elseif ((strncmp(s,'C',1)==1) & (length(s)==1)) & (((strncmp(u,u1,4)==1) & (bb1==bb) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'N',1)==1) & (length(s1)==1)) | ((strncmp(s1,'O',1)==1) & (length(s1)==1)) | ((strncmp(s1,'SG',2)==1) & (strncmp(bb,'CYS',3)==1)) | (strncmp(s1,'CG',2)==1) | (strncmp(s1,'CG1',3)==1) | (strncmp(s1,'CG2',3)==1) | ((strncmp(s1,'CD',2)==1) & (strncmp(bb,'PRO',3)==1)) | ((strncmp(s1,'OG',2)==1) & (strncmp(bb,'SER',3)==1)) | ((strncmp(s1,'OG1',3)==1) & (strncmp(bb,'THR',3)==1)))))
                    continue;
                elseif ((strncmp(s,'O',1)==1) & (length(s)==1)) & ((strncmp(u,u1,4)==1) & (bb1==bb) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'N',1)==1) & (length(s1)==1)) | ((strncmp(s1,'C',1)==1) & (length(s1)==1))))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'PRO',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CD',2)==1) & (strncmp(bb,'PRO',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'LYS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CE',2)==1) | (strncmp(s1,'NZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CD',2)==1) & (strncmp(bb,'LYS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CE',2)==1) | (strncmp(s1,'NZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CE',2)==1) & (strncmp(bb,'LYS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CD',2)==1) | (strncmp(s1,'CG',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'NZ',2)==1))
                    continue;
                elseif ((strncmp(s,'NZ',2)==1) & (strncmp(bb,'LYS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CD',2)==1) | (strncmp(s1,'CG',2)==1) | (strncmp(s1,'CE',2)==1))
                    continue;
                elseif ((strncmp(s,'CG1',3)==1) & (strncmp(bb,'ILE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CG2',3)==1))
                    continue;
                elseif ((strncmp(s,'CD1',3)==1) & (strncmp(bb,'ILE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG1',3)==1)  | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) |(strncmp(s1,'CG2',3)==1) )
                    continue;
                elseif ((strncmp(s,'CG2',3)==1) & (strncmp(bb,'ILE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG1',3)==1)  | (strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'N',1)==1) & (length(s1)==1)) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'PHE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD1',3)==1) |(strncmp(s1,'CD2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CD1',3)==1) & (strncmp(bb,'PHE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CD2',3)==1) & (strncmp(bb,'PHE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CE1',3)==1) & (strncmp(bb,'PHE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CE2',3)==1) & (strncmp(bb,'PHE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CZ',2)==1) & (strncmp(bb,'PHE',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CE1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD1',3)==1) |(strncmp(s1,'CD2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1) | (strncmp(s1,'OH',2)==1))
                    continue;
                elseif ((strncmp(s,'CD1',3)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1) | (strncmp(s1,'OH',2)==1) )
                    continue;
                elseif ((strncmp(s,'CD2',3)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1)  | (strncmp(s1,'OH',2)==1))
                    continue;
                elseif ((strncmp(s,'CE1',3)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ',2)==1)  | (strncmp(s1,'OH',2)==1))
                    continue;
                elseif ((strncmp(s,'CE2',3)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CZ',2)==1)  | (strncmp(s1,'OH',2)==1))
                    continue;
                elseif ((strncmp(s,'CZ',2)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'OH',2)==1))
                    continue;
                elseif ((strncmp(s,'OH',2)==1) & (strncmp(bb,'TYR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'CZ',2)==1) )
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD1',3)==1) |(strncmp(s1,'CD2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'NE1',3)==1))
                    continue;
                elseif ((strncmp(s,'CD1',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'NE1',3)==1))
                    continue;
                elseif ((strncmp(s,'CD2',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'NE1',3)==1))
                    continue;
                elseif ((strncmp(s,'NE1',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CE2',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'NE1',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CE3',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'NE1',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CZ2',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'NE1',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CZ3',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'NE1',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CH2',3)==1) | (strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CH2',3)==1) & (strncmp(bb,'TRP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1) | (strncmp(s1,'CE2',3)==1) | (strncmp(s1,'NE1',3)==1) | (strncmp(s1,'CE3',3)==1) | (strncmp(s1,'CZ2',3)==1) | (strncmp(s1,'CZ3',3)==1) | (strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'HIS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'ND1',3)==1) |(strncmp(s1,'CD2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'NE2',3)==1))
                    continue;
                elseif ((strncmp(s,'ND1',3)==1) & (strncmp(bb,'HIS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'NE2',3)==1))
                    continue;
                elseif ((strncmp(s,'CD2',3)==1) & (strncmp(bb,'HIS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'ND1',3)==1) | (strncmp(s1,'CE1',3)==1) | (strncmp(s1,'NE2',3)==1))
                    continue;
                elseif ((strncmp(s,'CE1',3)==1) & (strncmp(bb,'HIS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'ND1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'NE2',3)==1))
                    continue;
                elseif ((strncmp(s,'NE2',3)==1) & (strncmp(bb,'HIS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'ND1',3)==1) | (strncmp(s1,'CD2',3)==1) | (strncmp(s1,'CE1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'ASN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OD1',3)==1) |(strncmp(s1,'ND2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'OD1',3)==1) & (strncmp(bb,'ASN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'ND2',3)==1))
                    continue;
                elseif ((strncmp(s,'ND2',3)==1) & (strncmp(bb,'ASN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'OD1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'GLN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CD',2)==1) |  (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OE1',3)==1) |(strncmp(s1,'NE2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CD',2)==1) & (strncmp(bb,'GLN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |  (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OE1',3)==1) |(strncmp(s1,'NE2',3)==1))
                    continue;
                elseif ((strncmp(s,'OE1',3)==1) & (strncmp(bb,'GLN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) |(strncmp(s1,'NE2',3)==1))
                    continue;
                elseif ((strncmp(s,'NE2',3)==1) & (strncmp(bb,'GLN',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) |(strncmp(s1,'OE1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'GLU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CD',2)==1) |  (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OE1',3)==1) |(strncmp(s1,'OE2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CD',2)==1) & (strncmp(bb,'GLU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |  (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OE1',3)==1) |(strncmp(s1,'OE2',3)==1))
                    continue;
                elseif ((strncmp(s,'OE1',3)==1) & (strncmp(bb,'GLU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) |(strncmp(s1,'OE2',3)==1))
                    continue;
                elseif ((strncmp(s,'OE2',3)==1) & (strncmp(bb,'GLU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) |(strncmp(s1,'OE1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'ASP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OD1',3)==1) |(strncmp(s1,'OD2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'OD1',3)==1) & (strncmp(bb,'ASP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'OD2',3)==1))
                    continue;
                elseif ((strncmp(s,'OD2',3)==1) & (strncmp(bb,'ASP',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'OD1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'LEU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD1',3)==1) |(strncmp(s1,'CD2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CD1',3)==1) & (strncmp(bb,'LEU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD2',3)==1))
                    continue;
                elseif ((strncmp(s,'CD2',3)==1) & (strncmp(bb,'LEU',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD1',3)==1))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'MET',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'SD',2)==1) |(strncmp(s1,'CE',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'SD',2)==1) & (strncmp(bb,'MET',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CE',2)==1))
                    continue;
                elseif ((strncmp(s,'CE',2)==1) & (strncmp(bb,'MET',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'SD',2)==1))
                    continue;
                elseif ((strncmp(s,'OG',2)==1) & (strncmp(bb,'SER',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif (((strncmp(s,'SG',2)==1) & (strncmp(bb,'CYS',3)==1) & (strncmp(u,u1,3)~=1)) & (strncmp(s1,'SG',2)==1) & (strncmp(bb1,'CYS',3)==1) ) | ((strncmp(s,'SG',2)==1) & (strncmp(bb,'CYS',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'SG',2)==1))
                    continue;
                elseif ((strncmp(s,'OG1',3)==1) & (strncmp(bb,'THR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CG2',3)==1) & (strncmp(bb,'THR',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'OG1',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CG1',3)==1) & (strncmp(bb,'VAL',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG2',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CG2',3)==1) & (strncmp(bb,'VAL',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG1',3)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)))
                    continue;
                elseif ((strncmp(s,'CG',2)==1) & (strncmp(bb,'ARG',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & (((strncmp(s1,'N',1)==1) & (length(s1)==1)) | (strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CD',2)==1) |(strncmp(s1,'NE',2)==1) | ((strncmp(s1,'C',1)==1) & (length(s1)==1)) | (strncmp(s1,'CZ',2)==1))
                    continue;
                elseif ((strncmp(s,'CD',2)==1) & (strncmp(bb,'ARG',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CA',2)==1) | (strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'NE',2)==1) | (strncmp(s1,'CZ',2)==1) | (strncmp(s1,'NH1',3)==1) | (strncmp(s1,'NH2',3)==1))
                    continue;
                elseif ((strncmp(s,'NE',2)==1) & (strncmp(bb,'ARG',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CB',2)==1) | (strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD',2)==1) | (strncmp(s1,'CZ',2)==1) | (strncmp(s1,'NH1',3)==1) | (strncmp(s1,'NH2',3)==1))
                    continue;
                elseif ((strncmp(s,'CZ',2)==1) & (strncmp(bb,'ARG',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CG',2)==1) |(strncmp(s1,'CD',2)==1) | (strncmp(s1,'NE',2)==1) | (strncmp(s1,'NH1',3)==1) | (strncmp(s1,'NH2',3)==1))
                    continue;
                elseif ((strncmp(s,'NH1',3)==1) & (strncmp(bb,'ARG',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CD',2)==1) | (strncmp(s1,'NE',2)==1) | (strncmp(s1,'CZ',2)==1) | (strncmp(s1,'NH2',3)==1))
                    continue;
                elseif ((strncmp(s,'NH2',3)==1) & (strncmp(bb,'ARG',3)==1) & (strncmp(u,u1,4)==1) & (bb1==bb)) & ((strncmp(s1,'CD',2)==1) | (strncmp(s1,'NE',2)==1) | (strncmp(s1,'CZ',2)==1) | (strncmp(s1,'NH1',3)==1))
                    continue;
                elseif (strncmp(u,u1,4)==1) & (strncmp(bb,'HOH',3)==1)
                    continue;
                else
                    d2=dummy(g,1:3);
                    %aa2=pro_new(g,1:3);
                    if mean(d2)==0
                        continue;
                    end
                    h2=h_bond(g,1);
                    d3=sqrt((dummy(x,1)-dummy(g,1))^2 + (dummy(x,2)-dummy(g,2))^2 + (dummy(x,3)-dummy(g,3))^2);
                    if ((h1==1) && ((h2==3) || (h2==2))) || ((h1==2) && ((h2==1) || (h2==3))) || ((h1==3) && ((h2==1) || (h2==2))) || ((h1==3) && (h2==3))
                        if abs(d3) <= dist_h
                           % aa3p=sqrt((aa1p(1,1)-aa2(1,1))^2 + (aa1p(1,2)-aa2(1,2))^2 + (aa1p(1,3)-aa2(1,3))^2);
                            %if abs(d3) < aa3p
                                a=a+1;
                                %u_clash=u1;
                                %name_clash=bb1;
                                break;
                            %else continue;
                            %end
                        else continue;
                        end
                    elseif abs(d3) < dist_nb
                       % aa3p=sqrt((aa1p(1,1)-aa2(1,1))^2 + (aa1p(1,2)-aa2(1,2))^2 + (aa1p(1,3)-aa2(1,3))^2);
                       % if abs(d3) < aa3p
                            a=a+1;
                        %    u_clash=u1;
                         %   name_clash=bb1;
                            break;
                        %else continue;
                        %end
                    elseif ((strncmp(s,'SG',2)==1) && (strncmp(bb,'CYS',3)==1)) && ((strncmp(s1,'SG',2)==1) && (strncmp(bb1,'CYS',3)==1))
                        if (abs(aa3) >= 1.7) && (abs(aa3) <= 2.2)
                            continue;
                        end
                    end
                end
            end
        end
    end
    if a > 0
        break;
    end
end
out1 = a;
end
