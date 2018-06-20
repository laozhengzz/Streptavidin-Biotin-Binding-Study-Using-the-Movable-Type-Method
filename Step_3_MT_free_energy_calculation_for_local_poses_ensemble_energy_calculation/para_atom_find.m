function para_atom_final = para_atom_find(atom,mid_structure,input_structure)
                    
%atom = C2arXjz_atom;
%mid_structure = C2ar_atom;
%input_structure = C2arX_atom;
%oppo_atom1 = ON_atom;
%oppo_atom2 = ON_atom;


para_atom_final = [];





%%input_structure = protein;
torsion_mark = 0;
%%input_structure_Namelist = protein_Namelist;
%input_structure = protential_protein;
%input_structure(:,10) = 0;
tl=1;

candidate_atom=[];
candidate_atom = input_structure( dist(input_structure,atom)<=5 & dist(input_structure,atom)>0,:);



%a = 1;

%for j = 1:size(input_structure,1)
%    if norm(input_structure(j,1:3)-atom(1,1:3))<=5 && norm(input_structure(j,1:3)-atom(1,1:3))>0
%        candidate_atom(a,:) = input_structure(j,:);
        
%        a = a+1;
%    end
%end

mid_candidate_atom=[];
mid_candidate_atom = mid_structure( dist(mid_structure,atom)<=5 & dist(mid_structure,atom)>0,:);


%a = 1;

%for j = 1:size(mid_structure,1)
%    if sqrt((mid_structure(j,1)-atom(1,1))^2 + (mid_structure(j,2)-atom(1,2))^2 + (mid_structure(j,3)-atom(1,3))^2 )<=5 && sqrt((mid_structure(j,1)-atom(1,1))^2 + (mid_structure(j,2)-atom(1,2))^2 + (mid_structure(j,3)-atom(1,3))^2 )>0
%        mid_candidate_atom(a,:) = mid_structure(j,:);
%        a = a+1;
%    end
%end
des_atom_candi=[];
trig = 0;
for k = 1:size(mid_candidate_atom,1)
    if (norm(mid_candidate_atom(k,1:3) - atom(1,1:3))<2.0 && norm(mid_candidate_atom(k,1:3) - atom(1,1:3))>1 && (mid_candidate_atom(k,4)==20 || mid_candidate_atom(k,4)==21 || atom(1,4)==20 || atom(1,4)==21)) || (norm(mid_candidate_atom(k,1:3) - atom(1,1:3))<1.8 && norm(mid_candidate_atom(k,1:3) - atom(1,1:3))>1 && (mid_candidate_atom(k,4)~=20 && mid_candidate_atom(k,4)~=21 && atom(1,4)~=20 && atom(1,4)~=21))
        p_1 = mid_candidate_atom(k,:);
        %trig = 1;
        for l = 1:size(mid_candidate_atom,1)
            if ((norm(mid_candidate_atom(l,1:3) - p_1(1,1:3))<2.0 && norm(mid_candidate_atom(l,1:3) - p_1(1,1:3))>1 && (mid_candidate_atom(l,4)==20 || mid_candidate_atom(l,4)==21 || p_1(1,4)==20 || p_1(1,4)==21)) || (norm(mid_candidate_atom(l,1:3) - p_1(1,1:3))<1.8 && norm(mid_candidate_atom(l,1:3) - p_1(1,1:3))>1 && (mid_candidate_atom(l,4)~=20 && mid_candidate_atom(l,4)~=21 && p_1(1,4)~=20 && p_1(1,4)~=21))) && mid_candidate_atom(k,7)~=mid_candidate_atom(l,7)
                p_2 = mid_candidate_atom(l,:);
                %trig = 2;
                trig = 0;
                for m = 1:size(candidate_atom,1)
                    if ((norm(candidate_atom(m,1:3) - p_2(1,1:3))<2.0 && norm(candidate_atom(m,1:3) - p_2(1,1:3))>1 && (candidate_atom(m,4)==20 || candidate_atom(m,4)==21 || p_2(1,4)==20 || p_2(1,4)==21)) || (norm(candidate_atom(m,1:3) - p_2(1,1:3))<1.8 && norm(candidate_atom(m,1:3) - p_2(1,1:3))>1 && (candidate_atom(m,4)~=20 && candidate_atom(m,4)~=21 && p_2(1,4)~=20 && p_2(1,4)~=21))) &&  mid_candidate_atom(l,7)~=candidate_atom(m,7) && mid_candidate_atom(k,7)~=candidate_atom(m,7)   
                        p_3 = candidate_atom(m,:);
                        torsion_mark = torsion_mark+1;
                        des_atom_candi(tl,:)=p_3;
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



des_atom_candi = unique(des_atom_candi,'rows');

des_atom = des_atom_candi(dist(des_atom_candi,atom)<3.4 & dist(des_atom_candi,atom)>2.6,:);

para_atom_final = des_atom;

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






