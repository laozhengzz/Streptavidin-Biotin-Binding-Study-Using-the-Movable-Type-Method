function [Part_Matr_com,CNN1,VN1]  = intra_nonbond_potential_pp_GARFv1(protein,protein_Namelist,torsion_list_P,angle_list_P,nonb_d,PP_nonbonding_list,vdw_d,Rcutoff)

%torsion_list_P = torsion_list_P_all;
%angle_list_P = angle_list_P;
%nonb_d = vdw_d_all_AF;
%PP_nonbonding_list = vdw_all_AF_contact_list_num;
%vdw_d = nonb_d_f_refined;

protein_A = zeros(floor((1+size(protein,1))*size(protein,1)/2),10);
protein_B = zeros(floor((1+size(protein,1))*size(protein,1)/2),10);
protein_A_Namelist = zeros(floor((1+size(protein,1))*size(protein,1)/2),2);
protein_B_Namelist = zeros(floor((1+size(protein,1))*size(protein,1)/2),2);

pn=1;
for p1 = 1:size(protein,1)-1
    for p2 = p1+1:size(protein,1)
        if norm(protein(p1,1:3)-protein(p2,1:3))>=2.5 && norm(protein(p1,1:3)-protein(p2,1:3))<=Rcutoff
            protein_A(pn,:) = protein(p1,:);
            protein_B(pn,:) = protein(p2,:);
            protein_A_Namelist(pn,:) = protein_Namelist(p1,:);
            protein_B_Namelist(pn,:) = protein_Namelist(p2,:);
            
            pn=pn+1;
        end
    end
end

protein_A_Namelist(protein_A(:,7)==0,:) = [];
protein_B_Namelist(protein_B(:,7)==0,:) = [];

protein_A(protein_A(:,7)==0,:)=[];
protein_B(protein_B(:,7)==0,:)=[];


protein_A(:,11)=0;
protein_B(:,11)=0;

for i = 1:size(protein_A,1)
    for j = 1:size(torsion_list_P,1)
        if (protein_A(i,10) == torsion_list_P(j,10) && protein_B(i,10) == torsion_list_P(j,11)) || (protein_A(i,10) == torsion_list_P(j,11) && protein_B(i,10) == torsion_list_P(j,10))
            protein_A(i,11) = 1;
            protein_B(i,11) = 1;
            break
        end
    end
end


protein_A(protein_A(:,11)==1,:)=[];
protein_B(protein_B(:,11)==1,:)=[];


for i = 1:size(protein_A,1)
    for j = 1:size(angle_list_P,1)
        if (protein_A(i,7) == angle_list_P(j,8) && protein_B(i,7) == angle_list_P(j,9)) || (protein_A(i,7) == angle_list_P(j,9) && protein_B(i,7) == angle_list_P(j,8))
            protein_A(i,11) = 1;
            protein_B(i,11) = 1;
            break
        end
    end
end
protein_A_Namelist(protein_A(:,11)==1,:) = [];

protein_B_Namelist(protein_B(:,11)==1,:) = [];

protein_A(protein_A(:,11)==1,:)=[];
protein_B(protein_B(:,11)==1,:)=[];

protein_starterA = protein_A;
protein_starterB = protein_B;
protein_starterA_Namelist = protein_A_Namelist;
protein_starterB_Namelist = protein_B_Namelist;


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:size(protein_starterA,1)
    j = 0;
    lef=0;
    ind1 = [];
    ind2 = [];
    ind1 = find(PP_nonbonding_list(:,1) == protein_starterA_Namelist(i,1) & PP_nonbonding_list(:,2) == protein_starterA_Namelist(i,2) & PP_nonbonding_list(:,3) == protein_starterB_Namelist(i,1) & PP_nonbonding_list(:,4) == protein_starterB_Namelist(i,2));
    ind2 = find(PP_nonbonding_list(:,3) == protein_starterA_Namelist(i,1) & PP_nonbonding_list(:,4) == protein_starterA_Namelist(i,2) & PP_nonbonding_list(:,1) == protein_starterB_Namelist(i,1) & PP_nonbonding_list(:,2) == protein_starterB_Namelist(i,2));
    
    if ~isempty(ind1) || ~isempty(ind2)
        if isempty(ind1) && ~isempty(ind2)
            j = ind2;
        elseif ~isempty(ind1) && isempty(ind2)
            j = ind1;
        elseif ~isempty(ind1) && ~isempty(ind2)
            j = ind1;
        end
        
        dist=norm(protein_starterB(i,1:3)-protein_starterA(i,1:3));
        
        sd=find(abs(nonb_d(:,1)-dist)<=rd/2);
        if isempty(sd)
            continue
        end
        if (sd>halva) & (sd<size(nonb_d,1)-halva) %end
            com_nonb=nonb_d(sd-halva+1:sd+halva,j);
            %cf_vdw=repmat(vdw_f(sd-halva+1:sd+halva,order_num),round(800/(2*halva)),1);
        elseif sd>=size(nonb_d,1)-halva
            com_nonb=nonb_d(sd-(2*halva-1):sd,j);
            %cf_vdw=repmat(vdw_f(sd-(2*halva-1):sd,order_num),round(800/(2*halva)),1);
        elseif sd<=halva
            com_nonb=nonb_d(sd:sd+(2*halva-1),j);
            %cf_vdw=repmat(vdw_f(sd:sd+(2*halva-1),order_num),round(800/(2*halva)),1);
        end
        
        com_nonb_new=sum(com_nonb);
        
        if com_nonb_new>1
            Part_Matr_com=Part_Matr_com+log(com_nonb_new);
            lef=1;
        end
    end
        
    if (isempty(ind1) && isempty(ind2)) || lef==0
        
        if protein_starterA(i,4) >= protein_starterB(i,4)
            a = protein_starterB(i,4);
            b = protein_starterA(i,4);
        elseif protein_starterA(i,4) < protein_starterB(i,4)
            a = protein_starterA(i,4);
            b = protein_starterB(i,4);
        end
        
        dist=norm(protein_starterB(i,1:3)-protein_starterA(i,1:3)); %%cf_vdw_new=sum(cf_vdw);
        sd=find(abs(vdw_d(:,1)-dist)<=rd/2);   %%cf_vdw_new=1/cf_vdw_new;
        order_num=(22-(a-2) + 22)*(a-1)/2 + b-a+2;
        if length(sd)==0
            continue
        end
        if (sd>halva) & (sd<size(vdw_d,1)-halva) %end
            com_vdw=vdw_d(sd-halva+1:sd+halva,order_num);
            %cf_vdw=repmat(vdw_f(sd-halva+1:sd+halva,order_num),round(800/(2*halva)),1);
        elseif sd>=size(vdw_d,1)-halva
            com_vdw=vdw_d(sd-(2*halva-1):sd,order_num);
            %cf_vdw=repmat(vdw_f(sd-(2*halva-1):sd,order_num),round(800/(2*halva)),1);
        elseif sd<=halva
            com_vdw=vdw_d(sd:sd+(2*halva-1),order_num);
            %cf_vdw=repmat(vdw_f(sd:sd+(2*halva-1),order_num),round(800/(2*halva)),1);
        end
        
        
        com_vdw_new=sum(com_vdw);
        
        if com_vdw_new>0
            Part_Matr_com=Part_Matr_com+log(com_vdw_new);
            
        end
        
    else
        continue
    end
    
end


CNN1=size(protein_starterA,1);

VN1 = CNN1*log(const);





