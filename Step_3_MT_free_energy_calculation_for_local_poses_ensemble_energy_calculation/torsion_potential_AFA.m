%function [ZE_part_P,ZN_part_P] = torsion_potential(torsion_list,torsion_index,torsion)
function [Part_Matr_PT1,CNNP1,VL1,EPT1,  Part_Matr_PT2,CNNP2,VL2,EPT2] = torsion_potential_AFA(torsion_list1,torsion_list1_Namelist,torsion_list2,torsion_list2_Namelist,torsion_index,torsion,R_L_tor_srch)


%a=1;
%clear zn_atom;
%zn_atom=[];

%torsion_index(:,1:8) = floor(torsion_index(:,1:8));
%torsion_index2(:,1:4) = floor(torsion_index2(:,1:4));
%torsion_list1(:,1:8) = floor(torsion_list1(:,1:8));
%torsion_list2(:,1:8) = floor(torsion_list2(:,1:8));
%torsion_list1_1(:,1:4) = floor(torsion_list1_1(:,1:4));
%torsion_list2_2(:,1:4) = floor(torsion_list2_2(:,1:4));



Part_Matr_PT1=0;
Part_Matr_comA1=0;
Part_Matr_PT2=0;
Part_Matr_comA2=0;

EPT1=zeros(size(torsion_list1,1),1);
EijA1=zeros(size(torsion_list1,1),1);
EPT2=zeros(size(torsion_list2,1),1);
EijA2=zeros(size(torsion_list2,1),1);
%ligand initialization
%clear fidin;
%clear tline;

rd = 0.005;
hnstate = R_L_tor_srch/rd;
halva=hnstate;
halv=100;
halv1=halv/10;
const=halva*2;
const1=halv1*2;
column_num=100;

%for i = 2:size(torsion_index,1)
%    torsion_index_test1 = torsion_index{i,2};
%    torsion_index_test2 = torsion_index{i,4};
%    if strcmp(torsion_index_test1(1:2),'HI') == 1
%        torsion_index{i,2} = 'HI';
%    end
%    if strcmp(torsion_index_test2(1:2),'HI') == 1
%        torsion_index{i,4} = 'HI';
%    end
%end


pt =1 ;
for i = 1:size(torsion_list1,1)
    
    torsion_list1_Atom1 = torsion_list1_Namelist(i,1);
    torsion_list1_Res1 = torsion_list1_Namelist(i,2);
    torsion_list1_Atom2 = torsion_list1_Namelist(i,3);
    torsion_list1_Res2 = torsion_list1_Namelist(i,4);
    %if strcmp(torsion_list1_Res1(1:2),'HI') == 1
    %    torsion_list1_Res1 = 'HI';
    %end
    %if strcmp(torsion_list1_Res2(1:2),'HI') == 1
    %    torsion_list1_Res2 = 'HI';
    %end
    ind1 = [];
    ind2 = [];
    ind1 = find(torsion_index(:,1) == torsion_list1_Atom1 & torsion_index(:,2) == torsion_list1_Res1 & torsion_index(:,3) == torsion_list1_Atom2 & torsion_index(:,4) == torsion_list1_Res2);
    ind2 = find(torsion_index(:,3) == torsion_list1_Atom1 & torsion_index(:,4) == torsion_list1_Res1 & torsion_index(:,1) == torsion_list1_Atom2 & torsion_index(:,2) == torsion_list1_Res2);
    
    if ~isempty(ind1) || ~isempty(ind2)
        if isempty(ind1) && ~isempty(ind2)
            j = ind2;
        elseif ~isempty(ind1) && isempty(ind2)
            j = ind1;
        elseif ~isempty(ind1) && ~isempty(ind2)
            j = ind1;
        end
        sd = [];
        %torsion_list1(i,1:4) == torsion_index(j,1:4)% & ( (torsion_list(i,5) == torsion_list(i,6) & torsion_list(i,5) == torsion_list(i,7) & torsion_list(i,5) == torsion_list(i,8) & torsion_list(i,5) ~= 115 ) | (torsion_list(i,5) == torsion_list(i,6) & torsion_list(i,5) == torsion_list(i,7) & torsion_list(i,5) == torsion_list(i,8) & torsion_list(i,5) == 115 & torsion_index(j,5) == 115 ) )
        dist = torsion_list1(i,9);
        sd=find(abs(torsion(:,1)-dist)<=rd/2);
        
        if ~isempty(sd)
            
            if (sd>halva) && sd<length(torsion)-halva %end
                com_torsion=torsion(sd-halva+1:sd+halva,j);
                
            elseif sd>=length(torsion)-halva
                com_torsion=torsion(sd-(2*halva-1):sd,j);
                
            elseif sd<=halva
                com_torsion=torsion(sd:sd+(2*halva-1),j);
                
            end
            
            
            com_torsion_new=sum(com_torsion);
            
            
            if com_torsion_new>0
                %Part_Matr_PT=Part_Matr_PT+log(com_torsion_new);
                %com_torsion_new_A(i,1) = com_torsion_new;
                EPT1(i,1) = log(com_torsion_new);
                %pt=pt+1;
                
                
            end
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(torsion_list2,1)

    torsion_list2_Atom1 = torsion_list2_Namelist(i,1);
    torsion_list2_Res1 = torsion_list2_Namelist(i,2);
    torsion_list2_Atom2 = torsion_list2_Namelist(i,3);
    torsion_list2_Res2 = torsion_list2_Namelist(i,4);
    %if strcmp(torsion_list2_Res1(1:2),'HI') == 1
    %    torsion_list2_Res1 = 'HI';
    %end
    %if strcmp(torsion_list2_Res2(1:2),'HI') == 1
    %    torsion_list2_Res2 = 'HI';
    %end
    ind1 = [];
    ind2 = [];
    ind1 = find(torsion_index(:,1) == torsion_list2_Atom1 & torsion_index(:,2) == torsion_list2_Res1 & torsion_index(:,3) == torsion_list2_Atom2 & torsion_index(:,4) == torsion_list2_Res2);
    ind2 = find(torsion_index(:,3) == torsion_list2_Atom1 & torsion_index(:,4) == torsion_list2_Res1 & torsion_index(:,1) == torsion_list2_Atom2 & torsion_index(:,2) == torsion_list2_Res2);
    
    if ~isempty(ind1) || ~isempty(ind2)
        if isempty(ind1) && ~isempty(ind2)
            j = ind2;
        elseif ~isempty(ind1) && isempty(ind2)
            j = ind1;
        elseif ~isempty(ind1) && ~isempty(ind2)
            j = ind1;
        end
        sd = [];
        %torsion_list2(i,1:4) == torsion_index(j,1:4)% & ( (torsion_list(i,5) == torsion_list(i,6) & torsion_list(i,5) == torsion_list(i,7) & torsion_list(i,5) == torsion_list(i,8) & torsion_list(i,5) ~= 115 ) | (torsion_list(i,5) == torsion_list(i,6) & torsion_list(i,5) == torsion_list(i,7) & torsion_list(i,5) == torsion_list(i,8) & torsion_list(i,5) == 115 & torsion_index(j,5) == 115 ) )
        dist = torsion_list2(i,9);
        sd=find(abs(torsion(:,1)-dist)<=rd/2);
        
        if ~isempty(sd)
            
            if (sd>halva) && sd<length(torsion)-halva %end
                com_torsion=torsion(sd-halva+1:sd+halva,j);
                
            elseif sd>=length(torsion)-halva
                com_torsion=torsion(sd-(2*halva-1):sd,j);
                
            elseif sd<=halva
                com_torsion=torsion(sd:sd+(2*halva-1),j);
                
            end
            
            com_torsion_new=sum(com_torsion);
            
            
            if com_torsion_new>0
                %Part_Matr_PT=Part_Matr_PT+log(com_torsion_new);
                %com_torsion_new_A(i,1) = com_torsion_new;
                EPT2(i,1) = log(com_torsion_new);
                %pt=pt+1;
                
            end
        end
    end
    
end


for a = 1:size(EPT1,1)
    if EPT1(a,1)< -10 || EPT2(a,1)< -10 || EPT1(a,1)/EPT2(a,1) >100 || EPT2(a,1)/EPT1(a,1) >100
        EPT1(a,1)=0;
        EPT1(a,2)=-1;
        EPT2(a,1)=0;
        EPT2(a,2)=-1;
    end
end
%EPT1(EPT1(:,2)==-1,:)=[];
%EPT2(EPT2(:,2)==-1,:)=[];
Part_Matr_PT1 = sum(EPT1(:,1));
Part_Matr_PT2 = sum(EPT2(:,1));


if Part_Matr_PT1 ~= 0
    CNNP1=size(EPT1,1);
elseif Part_Matr_PT1 == 0
    CNNP1=0;
end
VL1=CNNP1*log(const);

if Part_Matr_PT2 ~= 0
    CNNP2=size(EPT2,1);
elseif Part_Matr_PT2 == 0
    CNNP2=0;
end
VL2=CNNP2*log(const);

