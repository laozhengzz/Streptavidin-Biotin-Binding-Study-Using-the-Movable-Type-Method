function result = check_planar(A,B,C,D,E)

%%%%%% check for planarity of atoms %%%%%%%%

x = (B-A);
y = (C-B);
z = (D-E);

norm11 = cross(x,y);
ans = dot(z,norm11);

if abs(ans) <= 0.20
    result = true;
else result = false;
end
    

end