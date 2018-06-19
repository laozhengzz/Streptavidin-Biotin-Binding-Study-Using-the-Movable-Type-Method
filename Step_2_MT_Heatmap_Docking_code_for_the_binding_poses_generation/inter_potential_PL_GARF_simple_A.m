function [Part_Matr_com,CNN1,VN1] = inter_potential_PL_GARF_simple_A(protein_starterA,protein_starterB,vdw_d)
rd = 0.005;
hnstate = 0.5/rd;
halva=hnstate;

const=halva*2;


Part_Matr_com=0;

for i=1:1:size(protein_starterA,1)
    
    a = protein_starterA(i,4);
    b = protein_starterB(i,4);
    
    dist=norm(protein_starterB(i,1:3)-protein_starterA(i,1:3)); 
    sd=[];
    sd=find(abs(vdw_d(:,1)-dist)<=rd/2);
    
    if isempty(sd)
        continue
    end
    
    order_num=(a-1)*34+b+1;
    if (sd>halva) & (sd<length(vdw_d)-halva)
        com_vdw=vdw_d(sd-halva+1:sd+halva,order_num);
        
    elseif sd>=length(vdw_d)-halva
        com_vdw=vdw_d(sd-(2*halva-1):sd,order_num);
        
    elseif sd<=halva
        com_vdw=vdw_d(sd:sd+(2*halva-1),order_num);
        
    end
    
    com_vdw_new=sum(com_vdw);
    
    if com_vdw_new>0
        Part_Matr_com=Part_Matr_com+log(com_vdw_new);
        
    end
    
    
    
end

CNN1=size(protein_starterA,1);

VN1 = CNN1*log(const);














