%%%%% make sure that the backbone is alternate %%%
%%%%% one angle is obtuse and two are acute    %%%
%%%%% check N,Ca and C atom position in backbone

function value = verifytriangle(A,B,C)

%%% measure all three angles and make sure that the triangle is obtuse
%%% with one angle greater than 90 degrees and two smaller than 90.
%%%        B
%%%        /\
%%%       /  \
%%%     A/____\C
%%%
%%% /_ABC > 90degrees (obtuse), /_BCA <= 50degree (acute), 
%%% /_BAC <= 50degrees (acute) 
    
%%% /_ABC
    x=(B-A);
    y=(C-B);

    theta=acos(dot(x,y)/(norm(x)*norm(y)));
    r1=radtodeg(theta);
    r11=180-r1;
    
 %%% /_ BCA
    x=(C-B); y = (A-C);
    theta=acos(dot(x,y)/(norm(x)*norm(y)));
    r2=rad2deg(theta);
    r22 = 180 -r2;
    
 %%% /_ BAC
    x=(A-B); y = (C-A);
    theta=acos(dot(x,y)/(norm(x)*norm(y)));
    r3=rad2deg(theta);
    r33 = 180 -r3;   
    
     if (r22 && r33) <= 37 && (r11 >= 108)
         value = 1;
     else value =0;
     end

end