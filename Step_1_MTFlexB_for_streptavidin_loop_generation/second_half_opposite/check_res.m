
 %%%%%%%%%%%%%%%%% checking vander waal collisions with rest of the protein %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
 
 
 function out = check_res(s,bb,u,dummy,k,i,name_atoms,hb)
 
if exist('chk_flip.mat')==2
     filename='chk_flip.mat';
      matfile(filename);
     
 end
 i_now=i; k_now=k; s_now=s; bb_now=bb; dummy_now=dummy; hb_now=hb; name_atoms_now=name_atoms;
load chk_flip.mat;
i=i_now; k=k_now; s=s_now; bb=bb_now; dummy=dummy_now; hb=hb_now; name_atoms=name_atoms_now;
clear i_now k_now s_now bb_now dummy_now hb_now name_atoms_now;


  dist_h=2.2; dist_nb=2.80;
aa=0;numn= nr(i,1); u_1=strtrim( adj_p(i,:)); u_11=strtrim( adj_n(i,:)); ngf=i;c_1=strtrim( adj_pc(i,:)); c_11=strtrim( adj_nc(i,:));
 
 aa1=dummy(k,1:3);
 mat_dist=pdist2( dummy,aa1);
 %aa1p= pro_new(l,1:3);
 h1= hb(k,1);
 
 %numnum=numn;
 asd= dummy;
 [k1,k2]=size(asd);
 if ((strncmp(s,'H1',2)==1) || (strncmp(s,'H2',2)==1)) && (strncmp(bb,'HOH',3)==1)
    dist_h = 1.4;
    dist_nb = 2.10;
end
 for ff=1:k1
     if (mat_dist(ff,1) < 2.8) && (ff ~= k)
         aa2= dummy(ff,1:3);
         s1=strtrim( name_atoms(ff,:));
         bb1=bb;
         u1=u;
         %c1=c;
         h2= hb(ff,1);
         if (strncmp(u_1,u1,4)==1) %& (strncmp(c_1,c1,3)==1)
             if (strncmp(bb1,'PRO',3)==1)
                 if ((strncmp(s,'N',1)==1) && (length(s)==1)) && ((strncmp(s1,'N',1)==1) && (length(s1)==1))
                     continue;
                 elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && (((strncmp(s1,'N',1)==1) && (length(s1)==1)) || (strncmp(s1,'CG',2)==1) || (strncmp(s1,'CA',2)==1) || (strncmp(s1,'CB',2)==1) || ((strncmp(s1,'C',1)==1) && (length(s1)==1)))
                     continue;
                 elseif ((strncmp(s,'C',1)==1) && (length(s)==1)) && (strncmp(s1,'CD',2)==1)
                     aa3=sqrt((aa1(1,1)-aa2(1,1))^2 + (aa1(1,2)-aa2(1,2))^2 + (aa1(1,3)-aa2(1,3))^2);
                     if aa3 >= 2.00
                         continue;
                     else aa=1; break;
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
             
         elseif (strncmp(u_11,u1,4)==1) %& (strncmp(c_11,c1,3)==1)
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
                else
                    
                    aa3=sqrt((aa1(1,1)-aa2(1,1))^2 + (aa1(1,2)-aa2(1,2))^2 + (aa1(1,3)-aa2(1,3))^2);
                    if ((h1==1) && ((h2==3) || (h2==2))) || ((h1==2) && ((h2==1) || (h2==3))) || ((h1==3) && ((h2==1) || (h2==2))) || ((h1==3) && (h2==3))
                        if abs(aa3) <= dist_h
                            %aa3p=sqrt((aa1p(1,1)-aa2(1,1))^2 + (aa1p(1,2)-aa2(1,2))^2 + (aa1p(1,3)-aa2(1,3))^2);
                           % if abs(aa3) < aa3p
                                aa=aa+1;
                                break;
                            %else continue;
                            %end
                        else continue;
                        end
                    elseif abs(aa3) < dist_nb
                        %aa3p=sqrt((aa1p(1,1)-aa2(1,1))^2 + (aa1p(1,2)-aa2(1,2))^2 + (aa1p(1,3)-aa2(1,3))^2);
                        %if abs(aa3) < aa3p
                            aa=aa+1;
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
 out = aa;
 clear k1 k2;
%  if aa > 0
%      clear temp; clear prob; clear energy;
%      temp(1:k,1:3)= pro_new(start:l,1:3); blb=1;
%      prob(1:k,1)=1; energy(1:k,1)=0;
%      dub(i,1)=1;n=2;
%      %break;
%  end
 
 clear mat_dist;
 end
