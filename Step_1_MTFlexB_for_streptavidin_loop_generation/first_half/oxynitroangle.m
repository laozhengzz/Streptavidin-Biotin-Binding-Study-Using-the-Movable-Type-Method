%%%%%% Make sure that backbone Oxygen is %%%%%%%%%%%%
%%%%%% at same angles from previous Ca and 
%%%%%% Nitrogen of next residue
%%%%%% The purpose is to make sure that 
%%%%%% Oxygen is in line with C



function value = oxynitroangle(A,B,C,D)

%%%%% /_ABC or /_Ca-C-O
    x=(B-A);
    y=(C-B);
    theta=acos(dot(x,y)/(norm(x)*norm(y)));
    r1=radtodeg(theta);
    r11=180-r1;

%%%%% /_DBC or /_N-C-O
    x=(B-D);
    y=(C-B);
    theta=acos(dot(x,y)/(norm(x)*norm(y)));
    r2=radtodeg(theta);
    r22=180-r2;

    
     if (r22 > 117 && r22 < 127) && (r11 > 117 && r11 < 127)
         value = 1;
     else value = 0;
     end

end