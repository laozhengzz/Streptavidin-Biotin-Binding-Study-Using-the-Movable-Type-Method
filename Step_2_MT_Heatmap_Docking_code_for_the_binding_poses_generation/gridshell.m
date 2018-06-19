function P_sh_grid = gridshell(lyr,centroid,protein_structure)

grid_d=0.4;
Rmax = lyr+2.5;
Rmax2 = lyr+2*grid_d;
max_x=0+Rmax2;
min_x=0-Rmax2;
max_y=0+Rmax2;
min_y=0-Rmax2;
max_z=0+Rmax2;
min_z=0-Rmax2;
sig_x=transpose(min_x:grid_d:max_x);
sig_y=transpose(min_y:grid_d:max_y);
sig_z=transpose(min_z:grid_d:max_z);

size_sig_x=size(sig_x);
size_sig_y=size(sig_y);
size_sig_z=size(sig_z);


for q=1:1:size_sig_x(1,1)
    grid_coor((q-1)*size_sig_y(1,1)*size_sig_z(1,1)+1:q*size_sig_y(1,1)*size_sig_z(1,1),1)=sig_x(q,1);
end
%toc
%tic
for q=1:1:size_sig_y(1,1)
    grid_coor1r((q-1)*size_sig_z(1,1)+1:q*size_sig_z(1,1),1)=sig_y(q,1);
end
%toc
%tic
grid_coor(:,2)=repmat(grid_coor1r,size_sig_x(1,1),1);
grid_coor(:,3)=repmat(sig_z,size_sig_x(1,1)*size_sig_y(1,1),1);
P_sh_num=0;
P_sh_grid_stand = zeros(10^5,3);
%for lyr = 2.7:grid_d:3.8
%for lyr = 2.7:grid_d:4    
    %clear P_sh
    P_sh=grid_coor(  (sqrt( (grid_coor(:,1)-0).^2+(grid_coor(:,2)-0).^2+(grid_coor(:,3)-0).^2 )< lyr+grid_d)&(sqrt( (grid_coor(:,1)-0).^2+(grid_coor(:,2)-0).^2+(grid_coor(:,3)-0).^2 )> lyr-grid_d/2),:);
    P_sh_grid_stand(P_sh_num+1:P_sh_num+size(P_sh,1),1:3)=P_sh;
    P_sh_num = P_sh_num+size(P_sh,1);
%end

P_sh_grid_stand(P_sh_grid_stand(:,1)==0,:)=[];
P_sh_grid_stand = unique(P_sh_grid_stand,'rows');

P_sh_grid(:,1) = P_sh_grid_stand(:,1) + centroid(1,1);
P_sh_grid(:,2) = P_sh_grid_stand(:,2) + centroid(1,2);
P_sh_grid(:,3) = P_sh_grid_stand(:,3) + centroid(1,3);

protein_structure( sqrt( (protein_structure(:,1)-centroid(1,1)).^2+(protein_structure(:,2)-centroid(1,2)).^2+(protein_structure(:,3)-centroid(1,3)).^2 )>Rmax,:)=[];


if size(protein_structure,1)>0
    for jz=1:1:size(protein_structure,1)
        P_sh_grid(sqrt( (P_sh_grid(:,1)-protein_structure(jz,1)).^2+(P_sh_grid(:,2)-protein_structure(jz,2)).^2+(P_sh_grid(:,3)-protein_structure(jz,3)).^2 )<=protein_structure(jz,6),:)=[];
    end
end











