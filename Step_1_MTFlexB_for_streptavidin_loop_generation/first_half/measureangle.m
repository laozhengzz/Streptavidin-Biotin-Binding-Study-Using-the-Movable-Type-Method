%%% measures angle b/w two vectors %%%

%%% It will be used to measure the angle between the growing loop by 
%%% calculating the angle between two known terminals : C of the i-1
%%% residue and N of the n+1 residue, where i-1 is the last residue before
%%% the loop starts and n+1 is the first residue after the loop ends. 
%%% The third point will be the C of the growing chain. So, basically the
%%% angle will be C(growing chain)-C (i-1)-N(n+1). If the angle is <= 100
%%% degrees, accept the move else reject it. 
%%% The purpose is to make sure that we are moving in the direction of the
%%% loop end.



function value = measureangle(A,B,C,count)
    x=(B-A);
    y=(C-B);

    theta=acos(dot(x,y)/(norm(x)*norm(y)));
    r1=radtodeg(theta);
    r2=180-r1;
     if (sign(count) == -1) && (r2 <= 90)  
         value = 1;
     elseif ((sign(count) == 1) || (sign(count) == 0)) && (r2 <= 130) %%% allow some slack for angle in the beginning of loop
         value = 1;
     else value = 0;
     end

end