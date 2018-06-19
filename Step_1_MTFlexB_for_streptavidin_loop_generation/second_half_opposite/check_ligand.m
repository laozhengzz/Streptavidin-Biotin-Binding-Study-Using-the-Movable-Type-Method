
 %%%%%%%%%%%%%%%%% checking vander waal collisions with rest of the protein %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
 
 
 function out3 = check_ligand(dummy,k,hbond)
 
 if exist('chk.mat')==2
     filename='chk.mat';
      matfile(filename);
     
 end
 k_now=k;  dummy_now=dummy; hbond_now=hbond;
load chk.mat;
k=k_now; dummy=dummy_now; hbond=hbond_now;
clear k_now dummy_now hbond_now;
  dist_h=2.1; dist_nb=2.8;
aa=0; 
 
 aa1=dummy(k,1:3);
 mat_dist=pdist2(lig_atom,aa1);
 h1= hbond;
 

 if ((strncmp(s,'H1',2)==1) || (strncmp(s,'H2',2)==1)) && (strncmp(bb,'HOH',3)==1)
    dist_h = 1.4;
    dist_nb = 2.10;
end
for ff=1:length(lig_atom)
    if mat_dist(ff,1) < 2.8
        aa2= lig_atom(ff,1:3);
        h2= lig_hbond(ff,1);
         
        aa3=sqrt((aa1(1,1)-aa2(1,1))^2 + (aa1(1,2)-aa2(1,2))^2 + (aa1(1,3)-aa2(1,3))^2);
        if ((h1==1) && ((h2==3) || (h2==2))) || ((h1==2) && ((h2==1) || (h2==3))) || ((h1==3) && ((h2==1) || (h2==2))) || ((h1==3) && (h2==3))
            if abs(aa3) <= dist_h
                aa=aa+1;
                break;
            else continue;
            end
        elseif abs(aa3) < dist_nb
            aa=aa+1;
            break;
            
        end
        
    end
end
out3 = aa;

 clear mat_dist;
 end
