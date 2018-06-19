function [repulse_num] = repulsenum(protein,ligand)
repulse_num=0;

for i=1:1:size(protein,1)
    for j=1:1:size(ligand,1)
        
        if (sqrt((ligand(j,1)-protein(i,1))^2+(ligand(j,2)-protein(i,2))^2+(ligand(j,3)-protein(i,3))^2)<2.5) && (((protein(i,4)==5 || protein(i,4)==11 || protein(i,4)==12 || protein(i,4)==13) && (ligand(j,4)>=7 && ligand(j,4)<=10)) || ((ligand(j,4)==5 || ligand(j,4)==11 || ligand(j,4)==12 || ligand(j,4)==13) && (protein(i,4)>=7 && protein(i,4)<=10)))                          
            repulse_num=repulse_num+1;
        elseif (sqrt((ligand(j,1)-protein(i,1))^2+(ligand(j,2)-protein(i,2))^2+(ligand(j,3)-protein(i,3))^2)<2.8)
            repulse_num=repulse_num+1;
            
        end
    end
end



