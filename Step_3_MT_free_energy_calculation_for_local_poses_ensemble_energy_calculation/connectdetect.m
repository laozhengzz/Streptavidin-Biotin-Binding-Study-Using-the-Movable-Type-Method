function connect_num = connectdetect(atom)

connect_num = 0;
for i = 1:size(atom,1)-1
    for j = i+1:size(atom,1)
        if norm(atom(i,1:3)-atom(j,1:3))<1.8
            connect_num = connect_num+1;
        end
    end
end