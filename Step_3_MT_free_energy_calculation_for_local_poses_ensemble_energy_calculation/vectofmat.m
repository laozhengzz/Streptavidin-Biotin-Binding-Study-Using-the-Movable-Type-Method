function logic = vectofmat(vec,mat)

logic = 0;
for i = 1:size(mat,1)
    if strfind(mat(i,:),vec)==1
        logic = 1;
        break
    end
end