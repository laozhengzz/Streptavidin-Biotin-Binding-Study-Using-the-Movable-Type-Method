%function [ZE_part_P,ZN_part_P] = torsion_potential(torsion_list,torsion_index,torsion)
function [Part_Matr_PT1,CNNP1,VL1,EPT1,Part_Matr_comA1,CNNA1,VNA1,EijA1] = torsion_potential_single(torsion_list1,torsion_list1_Namelist,torsion_index,torsion,R_L_tor_srch)


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

EPT1=zeros(size(torsion_list1,1),1);
EijA1=zeros(size(torsion_list1,1),1);

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
    ptn=0;
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
    
    
    for j = 1:size(torsion_index,1)
        
        if ( torsion_list1_Atom1==torsion_index(j,1) && torsion_list1_Atom2==torsion_index(j,3) && torsion_list1_Res1==torsion_index(j,2) && torsion_list1_Res2==torsion_index(j,4) ) || ( torsion_list1_Atom1==torsion_index(j,3) && torsion_list1_Atom2==torsion_index(j,1) && torsion_list1_Res1==torsion_index(j,4) && torsion_list1_Res2==torsion_index(j,2) )
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
                
                if length(sd)==0
                    break
                end
                com_torsion_new=sum(com_torsion);
                
                
                if com_torsion_new>0
                    %Part_Matr_PT=Part_Matr_PT+log(com_torsion_new);
                    %com_torsion_new_A(i,1) = com_torsion_new;
                    EPT1(i,1) = log(com_torsion_new);
                    %pt=pt+1;
                    ptn=1;
                    break
                end
            end
        end
        if ptn == 1
            break
            
        end
        
    end
    
end


%EPT1(EPT1(:,2)==-1,:)=[];
%EPT2(EPT2(:,2)==-1,:)=[];
Part_Matr_PT1 = sum(EPT1(:,1));


if Part_Matr_PT1 ~= 0
    CNNP1=size(EPT1,1);
elseif Part_Matr_PT1 == 0
    CNNP1=0;
end
VL1=CNNP1*log(const);






DOF_comA=0;
Part_Matr_comA=0;


for i = 1:size(torsion_list1,1)
    pt=0;
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
    
    
    for j = 1:size(torsion_index,1)
        if ( torsion_list1_Atom1==torsion_index(j,1) && torsion_list1_Atom2==torsion_index(j,3) && torsion_list1_Res1==torsion_index(j,2) && torsion_list1_Res2==torsion_index(j,4) ) || ( torsion_list1_Atom1==torsion_index(j,3) && torsion_list1_Atom2==torsion_index(j,1) && torsion_list1_Res1==torsion_index(j,4) && torsion_list1_Res2==torsion_index(j,2) )
            sd = [];
            dist = torsion_list1(i,9);
            sd=find(abs(torsion(:,1)-dist)<=rd/2);
            
            clear region_vdw
            if ~isempty(sd)
                
                if (sd>halva) & sd<length(torsion)-halva
                    region=torsion(sd-halva+1:sd+halva,j);
                elseif (sd<=halva)
                    region=torsion(sd:sd+(2*halva-1),j);
                elseif sd>=length(torsion)-halva
                    region=torsion(sd-(2*halva-1):sd,j);
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
                
                
                if length(sd)==0
                    continue
                end
                
                
                region_vdw_new=sum(region_vdw);
                
                
                region_vdw_new = sum(region_vdw_new')';
                
                
                if region_vdw_new>0
                    %Part_Matr_comA=Part_Matr_comA+log(region_vdw_new);
                    EijA1(i,1) = log(region_vdw_new);
                    pt=1;
                end
                if pt==1
                    break
                end
            end
        else
            continue;
        end
    end
    
    
end
%DOF_comA=sum(DOF_comA);
DOF_comA1 = size(EijA1,1);
CNNA1=DOF_comA1;

VNA1 = CNNA1*log(const1);


Part_Matr_comA1 = sum(EijA1(:,1));





