function [pocket] = pocket2find_PL_region(protein,ligand,Rcutoff)

m=1;
for i=1:1:size(protein,1)
    for j=1:1:size(ligand,1)
        if (sqrt((ligand(j,1)-protein(i,1))^2+(ligand(j,2)-protein(i,2))^2+(ligand(j,3)-protein(i,3))^2)<Rcutoff)
            protein_starterA(m,:)=protein(i,:);
            
            m=m+1;
            break
        end
    end
end

pocket = unique(protein_starterA,'rows');
pocket = sortrows(pocket,[7]); % rank according to the atomic index number



