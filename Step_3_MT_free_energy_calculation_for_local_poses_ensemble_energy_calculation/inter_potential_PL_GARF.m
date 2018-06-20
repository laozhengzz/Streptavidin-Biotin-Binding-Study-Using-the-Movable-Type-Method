%function [ZE_part_PL,ZE_part_PL2,ZN_part_PL,ZN_part_PL2] = inter_potential(protein_starterA,protein_starterB,vdw_d)
function [Part_Matr_com,Part_Matr_comA,CNN1,CNN2,VN1,VN2] = inter_potential_PL_GARF(protein_starterA,protein_starterB,vdw_d)
rd = 0.005;
hnstate = 0.2/rd;
halva=hnstate;
halv=100;
halv1=10;
const=halva*2;
const1=halv1*2;
column_num=100;
Part_Matr_com=0;
Part_Matr_comA=0;
for i=1:1:size(protein_starterA,1)
    
    a = protein_starterA(i,4);
    b = protein_starterB(i,4);
    
    dist=norm(protein_starterB(i,1:3)-protein_starterA(i,1:3)); %%cf_vdw_new=sum(cf_vdw);
    sd=[];
    sd=find(abs(vdw_d(:,1)-dist)<=rd/2);   %%cf_vdw_new=1/cf_vdw_new;
    
    if isempty(sd)
        continue
    end
    
    order_num=(a-1)*34+b+1;
    %c3n2
    if (sd>halva) & (sd<length(vdw_d)-halva) %end
        com_vdw=vdw_d(sd-halva+1:sd+halva,order_num);
        %cf_vdw=repmat(vdw_f(sd-halva+1:sd+halva,order_num),round(800/(2*halva)),1);
    elseif sd>=length(vdw_d)-halva
        com_vdw=vdw_d(sd-(2*halva-1):sd,order_num);
        %cf_vdw=repmat(vdw_f(sd-(2*halva-1):sd,order_num),round(800/(2*halva)),1);
    elseif sd<=halva
        com_vdw=vdw_d(sd:sd+(2*halva-1),order_num);
        %cf_vdw=repmat(vdw_f(sd:sd+(2*halva-1),order_num),round(800/(2*halva)),1);
    end
    
    com_vdw_new=sum(com_vdw);
    
    if com_vdw_new>0
        Part_Matr_com=Part_Matr_com+log(com_vdw_new);
        
        %Eij(i,1) = -0.5918*log(com_vdw_new);
    end
    
    
    
end




for i=1:1:size(protein_starterA,1)
    
    a = protein_starterA(i,4);
    b = protein_starterB(i,4);
    
    
    dist=norm(protein_starterB(i,1:3)-protein_starterA(i,1:3)); %%cf_vdw_new=sum(cf_vdw);
    
    sd=[];
    sd=find(abs(vdw_d(:,1)-dist)<=rd/2);   %%cf_vdw_new=1/cf_vdw_new;
    
    if isempty(sd)
        continue
    end
    order_num=(a-1)*34+b+1;
    %c3n2
    clear region_vdw
    
    if (sd>halva) && sd<length(vdw_d)-halva
        region=vdw_d(sd-halva+1:sd+halva,order_num);
    elseif (sd<=halva)
        region=vdw_d(sd:sd+(2*halva-1),order_num);
    elseif sd>=length(vdw_d)-halva
        region=vdw_d(sd-(2*halva-1):sd,order_num);
    end
    
    pr=0;
    
    size_region=size(region);
    for ri=6:size_region(1,1)-6
        if region(ri,1)>region(ri-1,1) && region(ri,1)>region(ri+1,1) && region(ri-1,1)>region(ri-2,1) && region(ri-2,1)>region(ri-3,1) && region(ri-3,1)>region(ri-4,1) && region(ri-4,1)>region(ri-5,1) && region(ri+1,1)>region(ri+2,1) && region(ri+2,1)>region(ri+3,1) && region(ri+3,1)>region(ri+4,1) && region(ri+4,1)>region(ri+5,1)
            pr=pr+1;
            if (ri>halv1) && ri<length(region)-halv1
                region_vdw(:,pr)=region(ri-halv1+1:ri+halv1,1);
            elseif (ri<=halv1)
                region_vdw(:,pr)=region(1:1+(2*halv1-1),1);
            elseif ri>=length(region)-halv1
                region_vdw(:,pr)=region(length(region)-(2*halv1-1):length(region),1);
            end
        end
    end
    if pr==0
        pr=pr+1;
        [ri,rii]=find(region==max(region));
        if (ri(1,1)>halv1) && ri(1,1)<length(region)-halv1
            region_vdw(:,pr)=region(ri(1,1)-halv1+1:ri(1,1)+halv1,1);
        elseif (ri(1,1)<=halv1)
            region_vdw(:,pr)=region(1:1+(2*halv1-1),1);
        elseif ri(1,1)>=length(region)-halv1
            region_vdw(:,pr)=region(length(region)-(2*halv1-1):length(region),1);
        end
    end
    
    DOF_comA(i,1)=pr;
    
    
    region_vdw_new=sum(region_vdw);
    
    
    
    if region_vdw_new>0
        Part_Matr_comA=Part_Matr_comA+log(region_vdw_new);
        %EijA(i,1) = -0.5918*log(region_vdw_new);
    end
    
    
end
DOF_comA=sum(DOF_comA);
Part_Matr_com = Part_Matr_com*1.4;
Part_Matr_comA = Part_Matr_comA*1.4;



CNN1=size(protein_starterA,1);
CNN2=DOF_comA;
ZN_part_PL=((sqrt(CNN1*8+1)+1)/2-4)*log(2)+((sqrt(CNN1*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN1*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN1*8+1)+1)/2-1)*log(const);
ZN_part_PL2=((sqrt(CNN2*8+1)+1)/2-4)*log(2)+((sqrt(CNN2*8+1)+1)/2-3)*log(2*pi)+((sqrt(CNN2*8+1)+1)/2-2)*log(4*pi)+((sqrt(CNN2*8+1)+1)/2-1)*log(const1);
            
ZE_part_PL=-0.5918*(Part_Matr_com-CNN1*log(const));
ZE_part_PL2=-0.5918*(Part_Matr_comA-DOF_comA*log(const));


VN1 = CNN1*log(const);
VN2 = CNN2*log(const1);














