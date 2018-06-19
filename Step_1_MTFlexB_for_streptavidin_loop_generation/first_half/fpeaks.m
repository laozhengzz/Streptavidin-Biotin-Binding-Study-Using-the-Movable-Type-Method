%%%%%%% find peaks in torsion values %%%%%%%%%%%%%%%%%%%%%%%
function [pro,d3,ener]=fpeaks(a,b,c)
%clear pro loc ene;
%%%%%% a ---> probability dataset, b----> energy dataset


[s1,s2]=size(a); x_prev=0;
p1=0;
for x=3:s1-2
    x1=a(x-2,1);
    x2=a(x-1,1);
    x3=a(x,1);
    x4=a(x+1,1);
    x5=a(x+2,1);

    if (x3 > x1) && (x3 > x2) && (x3 > x4) && (x3 > x5) && (x3 > 0.1) && (abs(x-x_prev) > 20) 
        p1=p1+1;
        pro(p1,1)=x3;
        loc(p1,1)=x;
        ener(p1,1)=b(x,1);
        d3(p1,1)=c(x,1);
        x_prev=x;   
    end
end
[m11,m12]=size(d3);
%%%%% blurring in peaks by selecting range of torsions %%%%%%%%%%%%
if (m11 <= 3) 
    clear x;
    aw=length(pro);
    ae=length(pro)+1;
    for x=1:aw
        loc(ae,1)=loc(x,1)-20;
        loc(ae+1,1)=loc(x,1)-10;
        loc(ae+2,1)=loc(x,1)+10;
        loc(ae+3,1)=loc(x,1)+20;
        pro(ae,1)=a(loc(ae,1),1);
        pro(ae+1,1)=a(loc(ae+1,1));
        pro(ae+2,1)=a(loc(ae+2,1),1);
        pro(ae+3,1)=a(loc(ae+3,1));
        ener(ae,1)=b(loc(ae,1),1);
        ener(ae+1,1)=b(loc(ae+1,1));
        ener(ae+2,1)=b(loc(ae+2,1));
        ener(ae+3,1)=b(loc(ae+3,1));
        d3(ae,1)=c(loc(ae,1),1);
        d3(ae+1,1)=c(loc(ae+1,1));
        d3(ae+2,1)=c(loc(ae+2,1));
        d3(ae+3,1)=c(loc(ae+3,1));
        ae=ae+4;
    end
end

