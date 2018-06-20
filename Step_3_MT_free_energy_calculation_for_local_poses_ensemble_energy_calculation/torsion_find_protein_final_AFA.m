function [candidate_Namelist_final,torsion_list] = torsion_find_protein_final_AFA(input_structure,input_structure_Namelist)
%%input_structure = protein;

%%input_structure_Namelist = protein_Namelist;
%input_structure = protential_protein;
%input_structure(:,10) = 0;
tl=1;

torsion_list = zeros(size(input_structure,1),11);
candidate_Namelist = zeros(size(input_structure,1),6);
for i = 1:size(input_structure,1)
    clear candidate_atom

    a = 1;
    candidate_atom(a,:) = input_structure(i,:);
    
    a = a+1;
    for j = 1:size(input_structure,1)
        if norm(input_structure(j,1:3)-input_structure(i,1:3))<=5 && i~=j
            candidate_atom(a,:) = input_structure(j,:);
            
            a = a+1;
        end
    end
    trig = 0;
    for k = 2:size(candidate_atom,1)
        if (norm(candidate_atom(k,1:3) - candidate_atom(1,1:3))<2.0 && (candidate_atom(k,4)==21 || candidate_atom(k,4)==22 || candidate_atom(1,4)==21 || candidate_atom(1,4)==22)) || (norm(candidate_atom(k,1:3) - candidate_atom(1,1:3))<1.8 && (candidate_atom(k,4)~=21 && candidate_atom(k,4)~=22 && candidate_atom(1,4)~=21 && candidate_atom(1,4)~=22))
            p_1 = candidate_atom(k,:);
            %trig = 1;
            for l = 2:size(candidate_atom,1)
                if ((norm(candidate_atom(l,1:3) - p_1(1,1:3))<2.0 && (candidate_atom(l,4)==21 || candidate_atom(l,4)==22 || p_1(1,4)==21 || p_1(1,4)==22)) || (norm(candidate_atom(l,1:3) - p_1(1,1:3))<1.8 && (candidate_atom(l,4)~=21 && candidate_atom(l,4)~=22 && p_1(1,4)~=21 & p_1(1,4)~=22))) && k ~=l
                    p_2 = candidate_atom(l,:);
                    %trig = 2;
                    trig = 0;
                    for m = 2:size(candidate_atom,1)
                        if ((norm(candidate_atom(m,1:3) - p_2(1,1:3))<2.0 && (candidate_atom(m,4)==21 || candidate_atom(m,4)==22 || p_2(1,4)==21 || p_2(1,4)==22)) || (norm(candidate_atom(m,1:3) - p_2(1,1:3))<1.8 && (candidate_atom(m,4)~=21 && candidate_atom(m,4)~=22 && p_2(1,4)~=21 & p_2(1,4)~=22))) && l ~=m && k ~=m
                            p_3 = candidate_atom(m,:);
                            candidate_Namelist(tl,1) = input_structure_Namelist(candidate_atom(1,10),1);
                            candidate_Namelist(tl,2) = input_structure_Namelist(candidate_atom(1,10),2);
                            candidate_Namelist(tl,3) = input_structure_Namelist(candidate_atom(m,10),1);
                            candidate_Namelist(tl,4) = input_structure_Namelist(candidate_atom(m,10),2);
                            candidate_Namelist(tl,5) = candidate_atom(1,10);
                            candidate_Namelist(tl,6) = p_3(1,10);
                            
                            %trig = 3;
                            %torsion_list(tl,1)=candidate_atom(1,7);
                            %torsion_list(tl,2)=p_1(1,7);
                            %torsion_list(tl,3)=p_2(1,7);
                            %torsion_list(tl,4)=p_3(1,7);
                            
                            torsion_list(tl,1)=candidate_atom(1,4);
                            torsion_list(tl,2)=p_1(1,4);
                            torsion_list(tl,3)=p_2(1,4);
                            torsion_list(tl,4)=p_3(1,4);
                            
                            torsion_list(tl,5)=candidate_atom(1,5);
                            torsion_list(tl,6)=p_1(1,5);
                            torsion_list(tl,7)=p_2(1,5);
                            torsion_list(tl,8)=p_3(1,5);
                            torsion_list(tl,9)=norm(candidate_atom(1,1:3) - p_3(1,1:3));
                            
                            torsion_list(tl,10)=candidate_atom(1,10);
                            torsion_list(tl,11)=p_3(1,10);
                            
                            
                            tl = tl+1;
                            trig = trig+1;
                            if trig ==3
                                break
                            end

                        end
                    end
                    %if trig ==1
                    %    break
                    %end


                end
            end
            %if trig ==1
            %    break
            %end

        end
    end


end

torsion_list(torsion_list(:,1)==0,:)=[];
candidate_Namelist(candidate_Namelist(:,1)==0,:)=[];
%id = cellfun('length',candidate_Namelist);
%candidate_Namelist(id(:,1)==0,:)=[];

for n = 1:size(torsion_list,1)
    %if torsion_list(n,1) > torsion_list(n,4) || (torsion_list(n,1) == torsion_list(n,4) && torsion_list(n,2) > torsion_list(n,3))
    if torsion_list(n,10) > torsion_list(n,11)
        temp_list = torsion_list(n,:);
        torsion_list(n,1) = temp_list(1,4);
        torsion_list(n,2) = temp_list(1,3);
        torsion_list(n,3) = temp_list(1,2);
        torsion_list(n,4) = temp_list(1,1);
        
        torsion_list(n,5) = temp_list(1,8);
        torsion_list(n,6) = temp_list(1,7);
        torsion_list(n,7) = temp_list(1,6);
        torsion_list(n,8) = temp_list(1,5);
        torsion_list(n,9) = temp_list(1,9);
        torsion_list(n,10) = temp_list(1,11);
        torsion_list(n,11) = temp_list(1,10);
        
        
        temp_Namelist = candidate_Namelist(n,:);
        candidate_Namelist(n,1) = temp_Namelist(1,3);
        candidate_Namelist(n,2) = temp_Namelist(1,4);
        candidate_Namelist(n,3) = temp_Namelist(1,1);
        candidate_Namelist(n,4) = temp_Namelist(1,2);
        candidate_Namelist(n,5) = temp_Namelist(1,6);
        candidate_Namelist(n,6) = temp_Namelist(1,5);
        
    end
end




torsion_list_index = unique(torsion_list(:,10:11),'rows');
[C,ia,ib] = intersect(torsion_list(:,10:11),torsion_list_index,'rows');
torsion_list = torsion_list(ia,:);

cnf = 1;

%candidate_Namelist_final = repmat({''},size(candidate_Namelist,1),4);
%candidate_Namelist_final(1:size(candidate_Namelist,1),5:6) = cell(size(candidate_Namelist,1),2);
%%%candidate_Namelist_final = cell(size(candidate_Namelist,1),6);
%%%for i = 1:size(torsion_list,1)
%%%    for j = 1:size(candidate_Namelist,1)
%%%        if candidate_Namelist{j,5} == torsion_list(i,10) && candidate_Namelist{j,6} == torsion_list(i,11)
%%%            candidate_Namelist_final(cnf,:) = candidate_Namelist(j,:);
%%%            cnf=cnf+1;
%%%            break
%%%        end
%%%    end
%%%end

[C,ia,ib] = intersect(torsion_list(:,10:11),candidate_Namelist(:,5:6),'rows');
candidate_Namelist_final = candidate_Namelist(ib,:);




%%%candidate_Namelist_final = cell(size(candidate_Namelist,1),6);
%%%for i = 1:size(torsion_list,1)
%%%    id = (cell2mat(candidate_Namelist(:,5)) == torsion_list(i,10) & cell2mat(candidate_Namelist(:,6)) == torsion_list(i,11));
%%%    candidate_Namelist_temp = candidate_Namelist(id,:);
%%%    candidate_Namelist_final(cnf,:) = candidate_Namelist_temp(1,:);
%%%    cnf=cnf+1;
%%%end
            
%%%id = cellfun('length',candidate_Namelist_final);
%%%candidate_Namelist_final(id(:,1)==0,:)=[];
%candidate_Namelist = uniqueRowsCA(candidate_Namelist);







