%%%%%% chirality determination %%%%%%%%%%%%%%%%%%%%%%

function out = chiral(chi_cen,ca,c3,c4)


%%%%%%%%%% translation of all the atoms %%%%%%%%%%%%%%%%
%%%    it will move chiral center to the origin      %%%

ca=ca-chi_cen;
c3=c3-chi_cen;
c4=c4-chi_cen;


%%%%%% make every atom as unit vector  %%%%%%%%%
ca=ca/norm(ca);
c3=c3/norm(c3);
c4=c4/norm(c4);

%%%%%%%%%% perform rotation on ca -i.e. 1st atom  %%%%%
%%%%% rotate other atoms using the same rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


GG=@(A,B)[dot(A,B),-norm(cross(A,B)),0;norm(cross(A,B)),dot(A,B),0;0,0,1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU = @(Fi,G) Fi*G*inv(Fi);

%%% rotating ca to the x-axis 
a = [-1 0 0]';   %%reference vector
b = ca';


U = UU(FFi(a,b), GG(a,b));  %%% rotation matrix

ca=ca*U;
c3=c3*U;
c4=c4*U;

%%%%%%% now rotating c3 keeping ca and chi_cen constant %%
a=[0.7071 0.7071 0]';
b=c3';
U = UU(FFi(a,b), GG(a,b)); 
c3=c3*U;
c4=c4*U;
out=sign(c4(1,3));
end