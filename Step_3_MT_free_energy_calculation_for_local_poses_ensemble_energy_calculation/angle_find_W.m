function angle_mark = angle_find(atom,mid_structure,input_structure)

%atom = ligand(jz,:);
%mid_structure = C2ar_atom;
%input_structure = Oall_atom;

angle_mark = 0;
%input_structure = protential_protein;
input_structure(:,10) = 0;
tl=1;
%for i = 1:size(input_structure,1)

candidate_atom=[];

a = 1;

for j = 1:size(input_structure,1)
    if sqrt((input_structure(j,1)-atom(1,1))^2 + (input_structure(j,2)-atom(1,2))^2 + (input_structure(j,3)-atom(1,3))^2 )<=4 && sqrt((input_structure(j,1)-atom(1,1))^2 + (input_structure(j,2)-atom(1,2))^2 + (input_structure(j,3)-atom(1,3))^2 )>0
        candidate_atom(a,:) = input_structure(j,:);
        a = a+1;
    end
end


mid_candidate_atom=[];

a = 1;

for j = 1:size(mid_structure,1)
    if sqrt((mid_structure(j,1)-atom(1,1))^2 + (mid_structure(j,2)-atom(1,2))^2 + (mid_structure(j,3)-atom(1,3))^2 )<=4 && sqrt((mid_structure(j,1)-atom(1,1))^2 + (mid_structure(j,2)-atom(1,2))^2 + (mid_structure(j,3)-atom(1,3))^2 )>0
        mid_candidate_atom(a,:) = mid_structure(j,:);
        a = a+1;
    end
end







trig = 0;
for k = 1:size(mid_candidate_atom,1)
    if sqrt((mid_candidate_atom(k,1)-atom(1,1))^2 + (mid_candidate_atom(k,2)-atom(1,2))^2 + (mid_candidate_atom(k,3)-atom(1,3))^2 )<2.1
        p_1 = mid_candidate_atom(k,:);
        trig = 0;
        for l = 1:size(candidate_atom,1)
            if sqrt((candidate_atom(l,1)-p_1(1,1))^2 + (candidate_atom(l,2)-p_1(1,2))^2 + (candidate_atom(l,3)-p_1(1,3))^2 )<2.1% && k ~=l
                p_2 = candidate_atom(l,:);
                %angle_list(tl,1)=candidate_atom(1,4);
                %angle_list(tl,2)=p_1(1,4);
                %angle_list(tl,3)=p_2(1,4);
                
                %angle_list(tl,4)=candidate_atom(1,5);
                %angle_list(tl,5)=p_1(1,5);
                %angle_list(tl,6)=p_2(1,5);
                %angle_list(tl,7)=norm(candidate_atom(1,1:3) - p_2(1,1:3));
                
                %angle_list(tl,8)=candidate_atom(1,7);
                %angle_list(tl,9)=p_2(1,7);
                
                angle_mark = angle_mark+1;
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


%end
















