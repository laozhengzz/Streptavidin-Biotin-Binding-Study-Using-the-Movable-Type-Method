
function mtflex_opposite_parallel(Begin,iRes)
Begin=str2num(Begin);

iRes=str2num(iRes);

example=matfile('tor.mat');
PR=example.P(:,Begin);
EN=example.E(:,Begin);
TO=example.Tor(:,Begin);

b=length(TO); a=1;
TO1=zeros(length(TO),1); EN1=TO1; PR1=TO1;

%%%%flip PR,EN and TO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:4
    d=3;
    aa1 = TO((b-d+1):b,:);
    aa2 = PR((b-d+1):b,:);
    aa3 = EN((b-d+1):b,:);
    a1=aa1(1,:);
    a2=aa1(3,:);
    a3=aa2(1,:);
    a4=aa2(3,:);
    a5=aa3(1,:);
    a6=aa3(3,:);
    aa1(1,:)=a2; aa1(3,:)=a1;
    aa2(1,:)=a4; aa2(3,:)=a3;
    aa3(1,:)=a6; aa3(3,:)=a5;
    TO1(a:(a+d-1),:)=aa1;
    EN1(a:(a+d-1),:)=aa3;
    PR1(a:(a+d-1),:)=aa2; 
 %   a=a+d; 
  %  b=b-d;
%end

clear TO EN PR;

TO=TO1; EN=EN1; PR=PR1;
clear TO1 EN1 PR1;

load chk.mat;
%%%%%%% Reverse everything %%%%%%%%%

%%%%% Point_back_n and Point_back_1
Point_back_11=Point_back_n;
Point_back_1n=Point_back_1;

Point_back_n=Point_back_1n;
Point_back_1=Point_back_11;

clear Point_back_11 Point_back_1n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%adj_p, adj_pc, adj_n, adj_nc, name_res
adj_p1=flipud(adj_p);
adj_n1=flipud(adj_n);
adj_nc1=flipud(adj_nc);
adj_pc1=flipud(adj_pc);
name_res1=flipud(name_res);
nr1=flipud(nr);
seq_res1=flipud(seq_res);
adj_p=adj_p1; adj_n=adj_n1;
adj_pc=adj_pc1; adj_nc=adj_nc1;
name_res=name_res1; nr=nr1;
seq_res=seq_res1;
clear adj_p1 adj_n1 adj_nc1 adj_pc1 name_res1 nr1 seq_res1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poc_coor1=zeros(size(poc_coor));

%%% Reverse the coordinates of the pocket %%%%%%%%%
a=1;
b=length(poc_coor);
for i=1:num_res
    d=nr(i,1);
    aa = poc_coor((b-d+1):b,:);
    a1=aa(1,:);
    a2=aa(3,:);
    aa(1,:)=a2; aa(3,:)=a1;
    clear a1 a2
    poc_coor1(a:(a+d-1),:)=aa;
    a=a+d; 
    b=b-d;
end

save chk_flip.mat;
original_pocket = name_res;

%%


%% declaring elements of protein with mutated side chains %%%
atom_new=char(zeros(size_pocket(1,1),3));
res_new=char(zeros(size_pocket(1,1),3));
seq_new=char(zeros(size_pocket(1,1),6));
chain  = char(zeros(size_pocket(1,1),1));
h_bond=zeros(size_pocket(1,1),1);
Num_pos_tot=zeros(length(nr),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MAXIMUM LENGTH OF THE LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is the maximum length possible for the entire chain
%%% It is decided based on the total number of residues plus
%%% the two additional N-C bonds with the known terminals.

Total_atoms  = num_res*3 + 2;

co=zeros(Total_atoms,1);
counter_atoms=zeros(Total_atoms,1);

c=1.340; b=1.540; a=1.480;

co(Total_atoms,1)=a;

for i = 2 : 3 : Total_atoms-1
    co(i,1) = a;
    co(i+1,1) = b;
    co(i+2,1) = c;
end

for i = 2 : Total_atoms
    counter_atoms(i,1)=sum(co(1:i,1));
end


match=zeros(num_res,1);
%%
%%%%%% Going over each residue one by one %%%%%%%%%%%%%%%%%%%%%%%%

bond_int = bonded_pro;
bond_energy = bonded_ene;
Mutated_pocket = original_pocket;  %%%defining the Mutated pocket, if needed


Final=0;   %% keeping tab of final number of PDBs
Col=1;
for r1=1:Col
    rn=1; l=0; 
    num_always=1; num1=1;
    drum11=[];
    l=0; y1=0; y2=1;
    Num_pos_tot=zeros(length(nr),3);
    st=1;
    count_back=1;
   fname=sprintf('chkpoint.mat');
    if exist(fname)==2
        load(fname);
        st=i+1;
        %%% Information needed from checkpoint file %%%%
        %%% We will only take iNum th conformation of the previously
        %%% generated conformations. So, we will keep only that
        %%% information.
         %d1=3*(iNum-1);
         %aaa=drum11(:,d1+1:d1+3);
         %t_n=temp_n(:,d1+1:d1+3); t_ca=temp_ca(:,d1+1:d1+3); t_c=temp_c(:,d1+1:d1+3);
         %t_o=temp_o(:,d1+1:d1+3);
         %PROB=pos_probability(:,iNum);
         %ENE=pos_energy(:,iNum);
         
         %clear pos_probabiity pos_energy drum11 temp_n temp_ca temp_c temp_o;
         %drum11=aaa; pos_probability = PROB; pos_energy=ENE; num1=1; num_always=1;
         %temp_n=t_n; temp_ca=t_ca; temp_c=t_c; temp_o=t_o;
         %clear aaa PROB ENE d1 t_n t_ca t_c t_o;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
i=iRes;
%i=iRes;
    bb=strtrim(Mutated_pocket(i,:));

    n=1; repeat = 0;
    
    if (strncmp(bb,bb1,3)==1)
        match(i,1)=0;   %% if User_seq and seq_pocket match, then no need of mutation   
        res_ori=poc_coor((sum(nr(1:i-1))+1):(sum(nr(1:i))),1:3);
    end
    res_ori=poc_coor((sum(nr(1:i-1))+1):(sum(nr(1:i))),1:3);
    [n_atoms, name_atoms, d1_atoms, d2_atoms, d3_atoms, hb] = aa_info_opp(bb);
    temp1=zeros(n_atoms,600000);

    pro1=ones(n_atoms,6000);
    ene1=zeros(n_atoms,6000);
    prob1=ones(n_atoms,6000);
    energy1=zeros(n_atoms,6000);
    clear tener tdum;
    num=num_always;
%     if (strncmp(bb,'GLY',3)==1) 
%         n_atoms=4;
%     else n_atoms=5;
%     end
    
    for k = 1:n_atoms
        clear A; clear B;
        clear Point1_1 Point2_1 Point3_1 d1 d2 d3;
        if k == 7
            p2=res_ori(2,1:3);
            p3=res_ori(5,1:3);
            p4=res_ori(6,1:3);
            p5=res_ori(7,1:3);
        end
        
        l=l+1;
        atom_new(l,:)=name_atoms(k,:);
        res_new(l,:)=Mutated_pocket(i,:);
        seq_new(l,:)=seq_res(i,:);
        chain(l,:)=chain_res(i,:);
        h_bond(l,1)=hb(k,1);
        
        s=strtrim(atom_new(l,:));
        u=strtrim(seq_new(l,:));
        
        if match(i,1) ~= 0
           temp1(k,1:3)=res_ori(k,1:3);
           prob1(k,1) = 1; energy1(k,1) = 0; 
        
        else  
        
            if (strncmp(s,'C',1)==1) && (length(s)==1)
                %a= bond_int(:,6);
                %b=bond_energy(:,6);
                %c=bond_int(:,1);
                %[pro,d3,ener]=fpeaks(a,b,c);
                pro=PR(rn,r1);
                d3=TO(rn,r1);
                ener=EN(rn,r1);
		count_back=count_back + 1;
                rn=rn+1; 
               
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                
                if i == 1
                    Point1_1=Point_back_1(1,:);
                    Point2_1=Point_back_1(2,:);
                    Point3_1=Point_back_1(3,:);
                    temp_o=Point_back_1(4,:);
                else
                    Point1_1=temp_n(1,:);
                    Point2_1=temp_ca(1,:);
                    
                    Point3_1=temp_c(1,:);    
                    
                end
                temp_try=temp_o;
                
            elseif (strncmp(s,'CA',2)==1)
                %a= bond_int(:,19);
                %b=bond_energy(:,19);
                %c=bond_int(:,1);
                %[pro,d3,ener]=fpeaks(a,b,c);
                pro=PR(rn,r1);
                d3=TO(rn,r1);
                ener=EN(rn,r1);
		count_back=count_back + 1;
                rn=rn+1;
                Point1_1=temp1(1,:);
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                temp_try=temp_o;
                
                if i == 1
                    [siz1,siz2]=size(Point1_1);
                    for h=1:(siz2/3)
                        Point2_1(1,(3*(h-1) + 1 : 3*(h-1) + 3))=Point_back_1(1,:);
                        Point3_1(1,(3*(h-1) + 1 : 3*(h-1) + 3))=Point_back_1(2,:);
                    end
                else
                    Point2_1=temp_n(1,:);
                    Point3_1=temp_ca(1,:);
                    
                end
                
                
            elseif (strncmp(s,'N',1)==1) && (length(s)==1)
                % a= bond_int(:,37);
                %b=bond_energy(:,37);
                %c=bond_int(:,1);
                %[pro,d3,ener]=fpeaks(a,b,c);
                pro=PR(rn,r1);
                d3=TO(rn,r1);
                ener=EN(rn,r1);
		count_back=count_back + 1;
                rn=rn+1;
                Point1_1=temp1(2,:);
                Point2_1=temp1(1,:);
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                temp_try=temp_o;
                if i == 1
                    [siz1,siz2]=size(Point1_1);
                    for h=1:(siz2/3)
                        Point3_1(1,(3*(h-1) + 1 : 3*(h-1) + 3))=Point_back_1(1,:);
                    end
                elseif i==num_res
                    [siz1,siz2]=size(Point1_1);
                    for h=1:(siz2/3)
                        Point3_1(1,(3*(h-1) + 1 : 3*(h-1) + 3))=Point_back_n(3,:);
                    end
                    d3=d1_atoms(1,1);
                    pro=1;
                    ener=0;
                else
                    Point3_1=temp_n(1,:);  
                end
                
                
            elseif ((strncmp(s,'O',1)==1) && (length(s)==1))         
                Point1_1=temp1(1,:);
                Point2_1=temp1(2,:);
                
                if i == 1
                    [siz1,siz2]=size(Point1_1);
                    for h=1:(siz2/3)
                        Point3_1(1,(3*(h-1) + 1 : 3*(h-1) + 3))=Point_back_1(1,:);
                    end
                else Point3_1=temp_n(1,:);
                    
                end
                d1=1.2336364;
                 d2=2.42;
%                 d3=2.28;

                d3=2.28; pro=1; ener=0;
                
                
                
            elseif (strncmp(s,'CB',2)==1)
                temp_try=temp1(4,:); %%%Oxygen
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(2,:);
                Point2_1=temp1(1,:);
                Point3_1=temp1(3,:);
                pro=1;
                ener=0;
                chi_cen=temp1(2,:); ca=temp1(3,:); c3=temp1(1,:);
                p2=res_ori(3,1:3);
                p3=res_ori(2,1:3);
                p4=res_ori(1,1:3);
                p5=res_ori(5,1:3);
                
                
            elseif ((strncmp(s,'CG',2)==1) && (length(s)==2)) || ((strncmp(s,'CG1',3)==1) && ((strncmp(bb,'ILE',3)==1) || (strncmp(bb,'VAL',3)==1)))
                
                if ((strncmp(bb,'PHE',3)==1) || (strncmp(bb,'TYR',3)==1))
                    a= bond_int(:,39);
                    b=bond_energy(:,39);
                    c=bond_int(:,1);
                    [pro,d3,ener]=fpeaks(a,b,c);
                    
                elseif ((strncmp(bb,'TRP',3)==1) || (strncmp(bb,'HIS',3)==1) || (strncmp(bb,'ASP',3)==1) || (strncmp(bb,'ASN',3)==1))
                    a= bond_int(:,38);
                    b=bond_energy(:,38);
                    c=bond_int(:,1);
                    [pro,d3,ener]=fpeaks(a,b,c);
                    
                elseif (strncmp(bb,'PRO',3)==1)
                    d3=d3_atoms(k,1);
                    pro=1;
                    ener=0;
                else
                    a= bond_int(:,30);
                    b=bond_energy(:,30);
                    c=bond_int(:,1);
                    [pro,d3,ener]=fpeaks(a,b,c);
                end
                
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1((k-1),:);
                Point2_1=temp1(2,:);
                Point3_1=temp1(1,:);
                
            elseif (strncmp(s,'CD',2)==1) && (strncmp(bb,'PRO',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                pro=1;
                ener=0;
                Point2_1=temp1((k-2),:);
                Point1_1=temp1((k-1),:);
                Point3_1=temp1(1,:);
                
            elseif (strncmp(s,'CG2',3)==1) && ((strncmp(bb,'ILE',3)==1) || (strncmp(bb,'VAL',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                pro=1;
                ener=0;
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                chi_cen=Point2_1; ca=Point3_1; c3=Point1_1;
                
            elseif (strncmp(s,'CD1',3)==1) && (strncmp(bb,'ILE',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                pro=1;
                ener=0;
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                
            elseif (strncmp(s,'CD2',3)==1) && (strncmp(bb,'LEU',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                pro=1;
                ener=0;
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                
            elseif ((strncmp(s,'CD',2)==1) && (strncmp(bb,'LYS',3)==1)) || ((strncmp(s,'CD1',3)==1) && (strncmp(bb,'LEU',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                a= bond_int(:,20);
                b=bond_energy(:,20);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'CE',2)==1) && (strncmp(bb,'LYS',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                a= bond_int(:,20);
                b=bond_energy(:,20);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'NZ',2)==1) && (strncmp(bb,'LYS',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                a= bond_int(:,26);
                b=bond_energy(:,26);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'SG',2)==1) && (strncmp(bb,'CYS',3)==1) && (bmb==0)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(2,:);
                Point3_1=temp1(3,:);
                a= bond_int(:,41);
                b=bond_energy(:,41);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
            elseif (strncmp(s,'ND1',3)==1) && (strncmp(bb,'HIS',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                pro=1;
                ener=0;
                
            elseif (strncmp(s,'CD2',3)==1) && (strncmp(bb,'HIS',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=Point3_1;
                B=Point2_1;
                C=Point1_1; E=B;
                
            elseif (strncmp(s,'CE1',3)==1) && (strncmp(bb,'HIS',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
                A=temp1(k-4,:);
                B=Point3_1;
                C=Point2_1;  E=C;
                
            elseif (strncmp(s,'NE2',3)==1) && (strncmp(bb,'HIS',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
                A=temp1(k-4,:);
                B=Point3_1;
                C=Point1_1; E=Point2_1;
                
            elseif (strncmp(s,'CD1',3)==1) && ((strncmp(bb,'TYR',3)==1)|| (strncmp(bb,'PHE',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                a= bond_int(:,17);
                b=bond_energy(:,17);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'CD2',3)==1) &&  ((strncmp(bb,'TYR',3)==1) || (strncmp(bb,'PHE',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
                A=Point3_1;
                B=Point1_1;
                C=Point2_1; E=Point1_1;
                
            elseif (strncmp(s,'CE1',3)==1) &&  ((strncmp(bb,'PHE',3)==1)|| (strncmp(bb,'TYR',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=temp1(k-4,:);
                B=Point3_1;
                C=Point1_1; E=Point1_1;
                
            elseif (strncmp(s,'CE2',3)==1) &&  ((strncmp(bb,'PHE',3)==1) || (strncmp(bb,'TYR',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=temp1(k-5,:);
                B=temp1(k-4,:);
                C=Point3_1; E=Point1_1;
                
            elseif (strncmp(s,'CZ',2)==1) &&  ((strncmp(bb,'PHE',3)==1) || (strncmp(bb,'TYR',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=temp1(k-6,:);
                B=Point3_1;
                C=temp1(k-4,:); E=Point1_1;
                
            elseif (strncmp(s,'OH',2)==1) &&  (strncmp(bb,'TYR',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=Point3_1;
                B=Point1_1;
                C=Point2_1; E=Point1_1;
                
            elseif (strncmp(s,'CD1',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
            elseif (strncmp(s,'CD2',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=Point3_1;
                B=Point2_1;
                C=Point1_1; E=Point1_1;
                
            elseif (strncmp(s,'NE1',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=temp1(k-4,:);
                B=Point3_1;
                C=Point2_1; E=Point2_1;
                
            elseif (strncmp(s,'CE2',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=temp1(k-5,:);
                B=temp1(k-4,:);
                C=Point3_1; E=Point2_1;
                
            elseif (strncmp(s,'CE3',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=Point3_1;
                B=Point1_1;
                C=Point2_1; E=Point3_1;
                
            elseif (strncmp(s,'CZ2',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=Point3_1;
                B=Point1_1;
                C=Point2_1; E=Point1_1;
                
            elseif (strncmp(s,'CZ3',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=Point1_1;
                B=Point3_1;
                C=Point2_1; E=Point1_1;
                
            elseif (strncmp(s,'CH2',3)==1) && (strncmp(bb,'TRP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                A=temp1(k-6,:);
                B=Point1_1;
                C=Point3_1; E=Point3_1;
                
            elseif (strncmp(s,'OD1',3)==1) && ((strncmp(bb,'ASN',3)==1)|| (strncmp(bb,'ASP',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                a= bond_int(:,29);
                b=bond_energy(:,29);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                
            elseif (strncmp(s,'ND2',3)==1) && (strncmp(bb,'ASN',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
            elseif (strncmp(s,'CD',2)==1) && ((strncmp(bb,'GLN',3)==1) || (strncmp(bb,'GLU',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                a= bond_int(:,10);
                b=bond_energy(:,10);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                
            elseif (strncmp(s,'OE1',3)==1) && ((strncmp(bb,'GLN',3)==1) || (strncmp(bb,'GLU',3)==1))
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                a= bond_int(:,29);
                b=bond_energy(:,29);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                
            elseif (strncmp(s,'NE2',3)==1) && (strncmp(bb,'GLN',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
            elseif (strncmp(s,'OE2',3)==1) && (strncmp(bb,'GLU',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
            elseif (strncmp(s,'OD2',3)==1) && (strncmp(bb,'ASP',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
                
            elseif (strncmp(s,'SD',2)==1) && (strncmp(bb,'MET',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                a= bond_int(:,43);
                b=bond_energy(:,43);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
            elseif (strncmp(s,'CE',2)==1) && (strncmp(bb,'MET',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                a= bond_int(:,21);
                b=bond_energy(:,21);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
            elseif (strncmp(s,'OG',2)==1) && (strncmp(bb,'SER',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(2,:);
                Point3_1=temp1(1,:);
                a= bond_int(:,40);
                b=bond_energy(:,40);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
            elseif (strncmp(s,'OG1',3)==1) && (strncmp(bb,'THR',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(2,:);
                Point3_1=temp1(1,:);
                a= bond_int(:,40);
                b=bond_energy(:,40);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'CG2',3)==1) && (strncmp(bb,'THR',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                pro=1;
                ener=0;
                chi_cen=Point2_1; ca=Point3_1; c3=Point1_1;
                
            elseif (strncmp(s,'CD',2)==1) && (strncmp(bb,'ARG',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(2,:);
                a= bond_int(:,10);
                b=bond_energy(:,10);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
            elseif (strncmp(s,'NE',2)==1) && (strncmp(bb,'ARG',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                a= bond_int(:,30);
                b=bond_energy(:,30);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
            elseif (strncmp(s,'CZ',2)==1) && (strncmp(bb,'ARG',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                a= bond_int(:,27);
                b=bond_energy(:,27);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'NH1',3)==1) && (strncmp(bb,'ARG',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                a= bond_int(:,24);
                b=bond_energy(:,24);
                c=bond_int(:,1);
                [pro,d3,ener]=fpeaks(a,b,c);
                
                
            elseif (strncmp(s,'NH2',3)==1) && (strncmp(bb,'ARG',3)==1)
                d1=d1_atoms(k,1);
                d2=d2_atoms(k,1);
                d3=d3_atoms(k,1);
                Point1_1=temp1(k-1,:);
                Point2_1=temp1(k-2,:);
                Point3_1=temp1(k-3,:);
                pro=1;
                ener=0;
            end
            
            
            D2=d2^2;
            
            n=1;
            
            [m5,m2]=size(d3); repeat=1;
            pro_orig=pro; d3_orig=d3; ener_orig=ener;
            
            while (n==1) && (repeat <=3)
                
                if (repeat == 1)
                    A_fpeaks=[d3_orig pro_orig ener_orig];
                    B_fpeaks=A_fpeaks;
                    [M,N]=size(B_fpeaks);
                    
                    
                    if (M > 5) && ((strncmp(s,'C',1)==1) && (length(s)==1)) %&& (i > 1)
                        d3=B_fpeaks(1:5,1); pro=B_fpeaks(1:5,2); ener=B_fpeaks(1:5,3);
                        d33=d3;
                        
                    elseif (M > 3) && (((strncmp(s,'N',1)==1) && (length(s)==1)) || (strncmp(s,'CA',2)==1))  %&& (i > 1)
                        d3=B_fpeaks(1:3,1); pro=B_fpeaks(1:3,2); ener=B_fpeaks(1:3,3);
                        d33=d3;
                        
                    elseif  (M > 3) && ((strncmp(s,'O',1)==1) && (length(s)==1))  %&& (i > 1)
                        d3=B_fpeaks(1,1); pro=B_fpeaks(1,2); ener=B_fpeaks(1,3);
                        d33=d3;
                        
                    elseif (M > 3) && (k > 5)  %&& (i > 1)
                        d3=B_fpeaks(1:3,1); pro=B_fpeaks(1:3,2); ener=B_fpeaks(1:3,3);
                        d33=d3;
                        
%                     elseif i == 1
%                         d3=B_fpeaks(:,1); pro=B_fpeaks(:,2); ener=B_fpeaks(:,3);
%                         d33=d3;
                        
                    else d3=B_fpeaks(:,1); pro=B_fpeaks(:,2); ener=B_fpeaks(:,3);
                        d33=d3;
                    end
                    
                    [m5,m2]=size(d33);
                    
                end
                
                if (repeat == 2)
                    if (M > 4)
                        d3=B_fpeaks(1:4,1); pro=B_fpeaks(1:4,2); ener=B_fpeaks(1:4,3);
                    else d3=B_fpeaks(:,1); pro=B_fpeaks(:,2); ener=B_fpeaks(:,3);
                    end
                    d33=d3;
                    [m5,m2]=size(d33);
                elseif (repeat >= 3)
                    d3=B_fpeaks(:,1); pro=B_fpeaks(:,2); ener=B_fpeaks(:,3);
                    d33=d3;
                    [m5,m2]=size(d33);
                end
                
		if (repeat > 1) && (M == 1) && (i <= (num_res-2)) && (k~=1)
                    d33=(linspace(d3-0.10,d3+0.10,3)'); pro=ones(3,1); ener=zeros(3,1);
                    [m5,m2]=size(d33);


            elseif (repeat == 2) && (M == 1) && (i <= (num_res-2)) && (k==1)
                   d33=[3.1850;3.4050;3.5300]; pro=[0.9643;0.5748; 0.5805]; ener=[-1.5571;-1.4040;1.4069];
                   [m5,m2]=size(d33);

                elseif (repeat > 2) && (M == 1) && (i <= (num_res-2)) && (k==1)
                d33=(linspace(d3-0.10,d3+0.10,3)'); pro=ones(3,1); ener=zeros(3,1);    
		[m5,m2]=size(d33);

                elseif (repeat > 1) && (M == 1) && ((i==num_res) || (i==num_res-1))
                    d33=(linspace(d3-0.20,d3+0.20,23)'); pro=ones(23,1); ener=zeros(23,1);
                    [m5,m2]=size(d33);
                end
               

                if repeat == 1
                    d31=d1; d32=d2;
                elseif repeat == 2
                    d31=(linspace(d1-0.01,d1+0.01,3)'); d32=d2;
                elseif repeat == 3
                    d31=(linspace(d1-0.02,d1+0.02,5)'); d32=(linspace(d2-0.01,d2+0.01,5)');
                elseif  repeat == 4
                    d31=(linspace(d1-0.01,d1+0.01,3)'); d32=(linspace(d2-0.01,d2+0.01,3)');
                elseif repeat == 5
                    d31=(linspace(d1-0.02,d1+0.02,5)'); d32=(linspace(d2-0.02,d2+0.02,5)');
                elseif repeat == 6
                    d31=(linspace(d1-0.03,d1+0.03,7)'); d32=(linspace(d2-0.03,d2+0.03,7)');
                elseif repeat == 7
                    d31=(linspace(d1-0.04,d1+0.04,9)'); d32=(linspace(d2-0.04,d2+0.04,9)');
                elseif repeat == 8
                    d31=(linspace(d1-0.05,d1+0.05,11)'); d32=(linspace(d2-0.05,d2+0.05,11)');
                elseif repeat == 9
                    d31=(linspace(d1-0.06,d1+0.06,13)'); d32=(linspace(d2-0.06,d2+0.06,13)');
                end
                
                for f=1:length(d31)
                    D1=d31(f)^2;
                    
                    for g=1:length(d32)
                        D2=d32(g)^2;
                        
                        for b=1:m5
                            D3=d33(b)^2;
                            prb=pro(b,1);
                            ene=ener(b,1);
                            d=0;
                            
                            for m=1:num
                                Point1=Point1_1(1,(d+1):(d+3));
                                Point2=Point2_1(1,(d+1):(d+3));
                                Point3=Point3_1(1,(d+1):(d+3));
                                if (k == 4)
                                    d1=1.2336364;
                                    alpha=2.4363;
                                    cos_beta=cos((2*pi -alpha)/2);
                                    a1=sqrt((Point2(1,1)-Point1(1,1))^2 + (Point2(1,2)-Point1(1,2))^2 + (Point2(1,3)-Point1(1,3))^2);
                                    a3=sqrt((Point3(1,1)-Point1(1,1))^2 + (Point3(1,2)-Point1(1,2))^2 + (Point3(1,3)-Point1(1,3))^2);
                                    d2=sqrt((a1)^2 + (d1)^2 -(2*a1*d1*cos_beta));
                                    d3=sqrt((a3)^2 + (d1)^2 -(2*a3*d1*cos_beta));  
                                end
                                
                                if (k == 1) || (k == 2) || (k == 3) %|| (k == 5)
                                    oo=temp_try(1,(d+1):(d+3));
                                elseif k == 5
                                    oo=temp1(4,(d+1):(d+3));
                                end
                                dummy=temp1((1:(k-1)),(d+1):(d+3));
                                if (strncmp(s,'O',1)==1) && (strncmp(bb,'HOH',3)==1)
                                    tdum=[];
                                    tener=[];
                                else
                                    tdum(1:k-1,1)=prob1(1:k-1,m);
                                    tener(1:k-1,1)=energy1(1:k-1,m);
                                end
                                [xf,yf,zf] = ana_solve(Point1,Point2,Point3,D1,D2,D3);
                                xf=real(xf); yf=real(yf); zf=real(zf);
                                
                                
                                for j = 1:2
                                    dummy(k,1:3) = [xf(j,1) yf(j,1) zf(j,1)];
                                    
                                    
                                    %%%%%% alternate backbone and check for position of Oxygen %%%%%%
                                    if ((strncmp(s,'N',1)==1) && (length(s)==1)) || ((strncmp(s,'C',1)==1) && (length(s)==1)) || (strncmp(s,'CA',2)==1) || ((strncmp(s,'O',1)==1) && (length(s)==1))
                                        value = verifytriangle(Point2,Point1,dummy(k,1:3));
                                    else value = 1;
                                    end
                                    %value = 1;
%                                     if ((strncmp(s,'N',1)==1) && (length(s)==1))
%                                         value2 = oxynitroangle(Point2,Point1,oo,dummy(k,1:3));
%                                     else value2 = 1;
%                                     end
                                    value2=1;
                                    
                                    value1=1;
                                    %%%%% planarity of rings %%%%
                                    if (exist('A')==1) && (exist('B')==1)
                                        result = check_planar(A(1,(d+1):(d+3)),B(1,(d+1):(d+3)),C(1,(d+1):(d+3)),dummy(k,1:3),E(1,(d+1):(d+3)));
                                    else result = true;
                                    end
                                    
                                    %%%% chirality of atom CA and CB %%%%
                                    if (strncmp(s,'CB',2)==1) || ((strncmp(s,'CG2',3)==1) && (strncmp(bb,'THR',3)==1)) || ((strncmp(s,'CG2',3)==1) && (strncmp(bb,'ILE',3)==1))
                                        pos_sign = chiral(chi_cen(1,(d+1):(d+3)),ca(1,(d+1):(d+3)),c3(1,(d+1):(d+3)),dummy(k,1:3));
                                        crystal_sign = chiral(p3,p2,p4,p5);
                                    else crystal_sign = 1; pos_sign = 1;
                                    end
                                    
                                    
                                    
                                    if n > 1
                                        for gi=1:n-1
                                            h2=3*(gi-1);
                                            f1=temp(:,h2+1:h2+3);
                                            f2=dummy;
                                            abb=((f1(:,1)-f2(:,1)).^2 + (f1(:,2)-f2(:,2)).^2 + (f1(:,3)-f2(:,3)).^2);
                                            euc_norm(1,gi)=sqrt(mean(abb));
                                        end
                                        %euc=(del(:,1)-dummy(k,1)).^2 + (del(:,2)-dummy(1,2)).^2 + (del(:,3)-dummy(1,3)).^2;
                                    else euc_norm = 1.0;
                                    end
                                    dis= sqrt((xf(j,1)-Point1(1,1))^2 + (yf(j,1)-Point1(1,2))^2 + (zf(j,1)-Point1(1,3))^2);
                                    dis1= sqrt((xf(j,1)-Point3(1,1))^2 + (yf(j,1)-Point3(1,2))^2 + (zf(j,1)-Point3(1,3))^2);
                                    dis2= sqrt((xf(j,1)-Point2(1,1))^2 + (yf(j,1)-Point2(1,2))^2 + (zf(j,1)-Point2(1,3))^2);
                                    
                                    if (strncmp(s,'CB',2)==1)
                                        dis3=sqrt((xf(j,1)-oo(1,1))^2 + (yf(j,1)-oo(1,2))^2 + (zf(j,1)-oo(1,3))^2);
                                    else dis3 = 2.00;
                                    end
                                    
%                                     if ((strncmp(s,'C',1)==1) && (length(s)==1))
%                                         dis4=sqrt((xf(j,1)-oo(1,1))^2 + (yf(j,1)-oo(1,2))^2 + (zf(j,1)-oo(1,3))^2);
%                                     else dis4 = 2.90;
%                                         
%                                     end
                                    dis4=2.90;
                                    if (strncmp(s,'CA',2)==1) && (i==num_res)
                                        dis5=sqrt((xf(j,1)-Point_back_n(1,1))^2 + (yf(j,1)-Point_back_n(1,2))^2 + (zf(j,1)-Point_back_n(1,3))^2);
                                    else dis5= 2.70;
                                    end
                                    
                                    
                                    %%% measure the length and compare to maximum length.
                                    if ((strncmp(s,'N',1)==1) && (length(s)==1)) || ((strncmp(s,'C',1)==1) && (length(s)==1)) || (strncmp(s,'CA',2)==1)
                                        
                                        dis6=sqrt((xf(j,1)-Point_back_n(1,1))^2 + (yf(j,1)-Point_back_n(1,2))^2 + (zf(j,1)-Point_back_n(1,3))^2);
                                        dis_t = counter_atoms(count_back,1);
                                        dis_remain = counter_atoms(Total_atoms,1) - dis_t;
                                        
                                        if dis6 <= dis_remain
                                            value3=1;
                                        else value3=0;
                                        end
                                        
                                    else value3=1;
                                    end
                                    
                                 %value3=1; 
                                 out = check_protein(s,bb,u,dummy,k,i,h_bond(l,1));
                                 out1 = check_pocket(s,bb,u,i,dummy,k,l,y1,h_bond(l,1));
                                 out2 = check_res(s,bb,u,dummy,k,i,name_atoms,hb);
                                 out3 = check_ligand(dummy,k,h_bond(l,1));
                                 %clash1 = check_back(h_back,coor_back,seq_back,res_back,atom_back,chain_back,adj_p,adj_n,adj_pc,adj_nc,i,dummy(k,1:3),s,bb,h_bond(l,1),u);
                              
                                 if (value == 1) && (value1 == 1) && (value2 == 1) && (value3 == 1)
                                     if (abs(dis) <= d31(f)+cut) && (abs(dis1) <= d33(b)+cut) && (abs(dis2) <= d32(g)+cut) && (abs(dis3) >= 2.0) && (abs(dis4) >= 2.8) && (abs(dis5) <= 2.70)
                                         if (out == 0) && (out1 == 0) && (out2 == 0) && (out3 == 0)
                                             if (crystal_sign == pos_sign) && result == true
                                                 if all(euc_norm > 0.01) == 1  
                                                       
                                                     if n == 1
                                                         temp=cat(2,dummy);
                                                         if k==1
                                                             temp_n=cat(2,Point1);
                                                             temp_ca=cat(2,Point2);
                                                             temp_o=cat(2,oo);
                                                         elseif k==2
                                                             %temp_ca=cat(2,xf(j,1));
                                                             temp_n=cat(2,Point2);
                                                             temp_o=cat(2,oo);
                                                          elseif k==3
                                                              temp_n=cat(2,Point3);
%                                                              temp_ca=cat(2,Point1);
%                                                              temp_c=cat(2,xf(j,1));
                                                         end
                                                         
                                                     else
                                                         temp=cat(2,temp,dummy);
                                                         if k==1
                                                             temp_n=cat(2,temp_n,Point1);
                                                             temp_ca=cat(2,temp_ca,Point2);
                                                             temp_o=cat(2,temp_o,oo);
                                                             
                                                         elseif k==2
                                                             %temp_ca=cat(2,temp_ca,xf(j,1));
                                                             temp_n=cat(2,temp_n,Point2);
                                                             temp_o=cat(2,temp_o,oo);
                                                          elseif k==3
                                                              temp_n=cat(2,temp_n,Point3);
%                                                              temp_ca=cat(2,temp_ca,Point1);
%                                                              temp_c=cat(2,temp_c,xf(j,1));
                                                         end
                                                         
                                                     end
                                                     tdum(k,1)=prb;
                                                     tener(k,1)=ene;
                                                     prob(1:k,n)=tdum; energy(1:k,n)=tener;
                                                     
                                                     n=n+1;
                                                 end
                                             end   
                                         end
                                     end
                                 end
                            end
                            d=3*m;
                        end
                    end
                end
            end
            repeat=repeat+1;
        end
        if (n == 1) && (repeat >= 4)
%             disp('ERROR: This residue cannot fit in the binding pocket:');
%             i, bb
            break;
        end
        
        temp1=temp;
        num=n-1; prob1=prob; energy1=energy;
        clear temp prob energy;
        
        end
        
    end
    if (n == 1) && (repeat >= 4)
        break;
        %disp('Please choose another residue at this position. Quit code now by pressing ctrl+c');
        %disp('<a href="MATLAB: dbquit;">Yes</a> / <a href="MATLAB: dbcont;">No</a>');
        keyboard;
        %break;
    end
    
    temp1=temp1(:,1:3*num);
    prob1=prob1(:,1:num);
    energy1=energy1(:,1:num);
    
    pos_prob = (prod(prob1,1))';    %% probability of each pose for ith residue
    pos_ener = (sum(energy1,1))';   %% energy of each pose for ith residue
    
    Num_pos_tot(i,1) = num;
    
    
    
    %%%%%%%% Keeping poses certain distance apart and killing rest %%%%%%%%
    if (num > 1e4) && (num <= 1e5)
        dem=0.3;
    elseif (num > 1e3) && (num <= 1e4)
        dem=0.2;
    elseif (num > 1e2) && (num <= 1e3)
        dem=0.1;
    elseif (num > 1e5)
        dem=0.5;
    elseif (num <= 1e2)
        dem=0.05;
    else dem=0.8;
    end
    
%     if (i == 1) || (i == num_res -1) || (i == num_res)
% 	dem=0.3;
% 	else dem=0.5;
%     end
    dem=0.20;
	if dem ~= 0
        o1=0;
        
        if num > 1
            clear pos_prob1 pos_ener1 temp11;
            pos_prob1=ones(size(pos_prob)); pos_ener1=zeros(size(pos_prob));
            siz=1; rr(1,1)=1;
            o1=1; o2=3*(o1-1);
            siz=o1; rr(o1,1)=1;
            temp11(:,o2+1:o2+3)=temp1(:,1:3);
            pos_prob1(o1,1)=pos_prob(1,1);
            pos_ener1(o1,1)=pos_ener(1,1);
            for gi=2:num
                h2=3*(gi-1);
                f1=temp1(:,h2+1:h2+3);
                ab=ones(1,siz);
                for gg=1:siz
                    h1=3*(rr(gg,1)-1);
                    f2=temp1(:,h1+1:h1+3);
                    abb=((f1(:,1)-f2(:,1)).^2 + (f1(:,2)-f2(:,2)).^2 + (f1(:,3)-f2(:,3)).^2);
                    ab(1,gg)=sqrt(mean(abb));
                end
                if  (any(ab < dem) == 0) 
                    o1=o1+1; o2=3*(o1-1);
                    siz=o1; rr(o1,1)=gi;
                    temp11(:,o2+1:o2+3)=f1;
                    pos_prob1(o1,1)=pos_prob(gi,1);
                    pos_ener1(o1,1)=pos_ener(gi,1);
                end
                
            end
            
            pos_prob1=pos_prob1(1:o1,1); pos_ener1=pos_ener1(1:o1,1);
            
            num=o1;
            clear temp1 pos_prob pos_ener;
            %pos_prob = PROB1; pos_ener = ENE1;
            temp1 = temp11; pos_prob=pos_prob1; pos_ener=pos_ener1;
            clear temp11 pos_prob1 pos_ener1;
        end
        
    end
    [size_drum11,f]=size(drum11); %%% get the size of drum11 matrix
    
    %%%%%%%%%%%%%%%%  permuting and combining poses of different residues  %%%%%%%%%%%%%%%%
   clear i1;
    
    s1=strtrim(atom_new(y1+1,:));
    u1=strtrim(seq_new(y1+1,:));
    c1=strtrim(chain(y1+1,:));
    u_1=strtrim(adj_p(i,:));
    if i > 1
        for j=l:-1:1
            s=strtrim(atom_new(j,:)); u=strtrim(seq_new(j,:));
            if ((strncmp(s,'N',1)==1) && (length(s)==1)) && (strncmp(u_1,u,4)==1)
                i1=j;
                break;
            end
        end
    end
    
%     if i >= 4
%         u_1=strtrim(adj_n(i-2,:));
%         for j=l:-1:1
%             s=strtrim(atom_new(j,:)); u=strtrim(seq_new(j,:));
%             if (strncmp(s,'CA',2)==1) && (strncmp(u_1,u,4)==1)
%                 pos_ca=j;
%                 break;
%             end
%         end
%         
%     end
    
    
    j2 = 0;
    repeats = 0;
%     clear temp_n temp_c temp_ca temp_o;
%     clear diff_n1 doff_o1 diff_ca1 diff_c1 diff;
%     temp_n1=[0 0 0]; temp_c1=temp_n1; temp_ca1=temp_n1; temp_o1=temp_n1;
    count=0;
    for f=1:num1
        d1=3*(f-1);
        if i ~= 1
            aaa=drum11(:,d1+1:d1+3);
            PROB=pos_probability(1:i-1,f);
            ENE=pos_energy(1:i-1,f);
        else
            aaa=[]; PROB=[]; ENE=[];
        end
        for m=1:num
            d=3*(m-1);
            e1=temp1(:,d+1:d+3);
            if (i ~= 1) && (exist('i1','var') ~= 0)
                i11=e1(1,:); i12=aaa(i1,:); i13=aaa(i1-1,:);i14=e1(2,:); i15=e1(4,:);
                aa3p=sqrt((i11(1,1)-i12(1,1))^2 + (i11(1,2)-i12(1,2))^2 + (i11(1,3)-i12(1,3))^2);
                aa3pp=sqrt((i11(1,1)-i13(1,1))^2 + (i11(1,2)-i13(1,2))^2 + (i11(1,3)-i13(1,3))^2);
                aa3ppp=sqrt((i12(1,1)-i14(1,1))^2 + (i12(1,2)-i14(1,2))^2 + (i12(1,3)-i14(1,3))^2);
                aa3pppp=sqrt((i12(1,1)-i15(1,1))^2 + (i12(1,2)-i15(1,2))^2 + (i12(1,3)-i15(1,3))^2); 
                
            else aa3p=1.5; aa3pp=2.0; aa3ppp=2.0; aa3pppp=2.0;
            end
%             if (i >= 4) && (exist('pos_ca','var') ~= 0)
%                 curr_ca=e1(2,:); old_ca=aaa(pos_ca,:);
%                 meas_cadist=sqrt((curr_ca(1,1)-old_ca(1,1))^2 + (curr_ca(1,2)-old_ca(1,2))^2 + (curr_ca(1,3)-old_ca(1,3))^2);
%else
    meas_cadist=13.99;
 %           end
            if (aa3p < 1.7) && (aa3pp >= 2.0) && (aa3ppp >= 2.0) && (aa3pppp >= 2.0) && (meas_cadist <= 14.00)
                drume=[aaa;e1];
                out1 = check_pocket_v2(h_bond,atom_new,res_new,seq_new,chain,i,drume,l,y1);
                
                if out1 == 0
%                     %%% to find the unique N-Ca-C-O combinations only
%                     diff_n1=pdist2(temp_n1,e1(1,:));
%                     diff_ca1=pdist2(temp_ca1,e1(2,:));
%                     diff_c1=pdist2(temp_c1,e1(3,:));
%                     diff_o1=pdist2(temp_o1,e1(4,:));
%                     diff= diff_n1 + diff_ca1 + diff_c1 + diff_o1;
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    j1=3*j2;
                    drum(:,j1+1:j1+3)=drume;     
                    
                    j2=j2+1;
                    PROB1(:,j2)=[PROB;pos_prob(m,1)];
                    ENE1(:,j2)=[ENE;pos_ener(m,1)];
                    
%                     %if (any(diff (:) == 0) == 0)
%                         c=3*count;
%                         temp_n(:,c+1:c+3)=e1(1,:);
%                         temp_ca(:,c+1:c+3)=e1(2,:);
%                         temp_c(:,c+1:c+3)=e1(3,:);
%                         temp_o(:,c+1:c+3)=e1(4,:);
%                         
%                         count=count+1;
%                         %%to keep tab of uncommon backbones
%                         temp_n1(count,:)=e1(1,:);
%                         temp_ca1(count,:)=e1(2,:);
%                         temp_c1(count,:)=e1(3,:);
%                         temp_o1(count,:)=e1(4,:);
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                     %end
                end
            end
        end
    end
   
    num1 = j2;
    
    if num1 == 0
        break;
    else
    clear drume PROB pos_prob ENE pos_ener temp_n1 temp_ca1 temp_c1 temp_o1;
    pos_prob = PROB1; pos_ener = ENE1; drume = drum;
    %clear PROB1 ENE1;
    
    Num_pos_tot(i,2) = num1;
    
    num_always=count;  %%unique combinations of backbone
    end
    if i > 1 
        sss=sum(pos_ener);
    
        sss1=[pos_ener;sss];
    else sss1=pos_ener; sss=sss1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delR1=sort(sss,'ascend');
clear drum PROB1 ENE1;
xx=ind2sub(length(pos_prob),find(sss1(end,:)==delR1(1,1)));

ay=0; t=1;
while t <= num1
     aw=0;
     xx=ind2sub(length(pos_prob),find(sss1(end,:)==delR1(1,t)));
     %xx=ind2sub(length(pos_prob),find(pos_prob(1,:)==delR1(1,t)));
     [xx1,xx2]=size(xx);
    
    while aw < xx2
        aw=aw+1;
        ay=ay+1;
        j=3*(ay-1); tr=3*(xx(1,aw)-1);
        PROB1(:,ay)=pos_prob(:,xx(1,aw));
        ENE1(:,ay) = pos_ener(:,xx(1,aw));
        drum(:,j+1:j+3)=drume(:,tr+1:tr+3);
        t=t+1;
%         if aw==xx2
%             aw=0;
%         end
    end

 end

num1=ay;
clear pos_prob pos_ener drume;
pos_prob=PROB1; pos_ener=ENE1;



dem=0.20;

if (dem ~= 0) && (num1 > 1)
    o1=0;
    clear pos_prob pos_ener drume temp11;
    clear temp_n temp_c temp_ca temp_o;
    clear diff_n1 doff_o1 diff_ca1 diff_c1 diff;
    temp_n1=[0 0 0]; temp_c1=temp_n1; temp_ca1=temp_n1; temp_o1=temp_n1;
    count=0;
    
    %if num1 > 1
    clear pos_prob pos_ener;
    o1=1; o2=3*(o1-1);
    siz=o1; rr(o1,1)=1;
    temp11(:,o2+1:o2+3)=drum(1:l,1:3);
    pos_prob(:,o1)=PROB1(:,1);
    pos_ener(:,o1)=ENE1(:,1);
    temp_c1=temp11(size_drum11+1,1:3);
    temp_ca1=temp11(size_drum11+2,1:3);
    temp_n1=temp11(size_drum11+3,1:3);
    temp_o1=temp11(size_drum11+4,1:3);
    temp_n=temp_n1; temp_ca=temp_ca1; temp_c=temp_c1; temp_o=temp_o1;
    count=1;
    for gi=2:num1
        h2=3*(gi-1);
        f1=drum(1:l,h2+1:h2+3);
        ab=ones(1,siz);
        for gg=1:siz
            h1=3*(rr(gg,1)-1);
            f2=drum(1:l,h1+1:h1+3);
            abb=((f1(:,1)-f2(:,1)).^2 + (f1(:,2)-f2(:,2)).^2 + (f1(:,3)-f2(:,3)).^2);
            ab(1,gg)=sqrt(mean(abb));
        end
        if (any(ab < dem) == 0)
            o1=o1+1; o2=3*(o1-1);
            siz=o1; rr(o1,1)=gi;
            temp11(:,o2+1:o2+3)=f1;
            pos_prob(:,o1)=PROB1(:,gi);
            pos_ener(:,o1)=ENE1(:,gi);
            
            %%% to find the unique N-Ca-C-O combinations only
            diff_c1=pdist2(temp_c1,f1(size_drum11+1,:));
            diff_ca1=pdist2(temp_ca1,f1(size_drum11+2,:));
            diff_n1=pdist2(temp_n1,f1(size_drum11+3,:));
            diff_o1=pdist2(temp_o1,f1(size_drum11+4,:));
            diff= diff_n1 + diff_ca1 + diff_c1 + diff_o1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if (any(diff (:) == 0) == 0)
                c=3*count;
                temp_c(:,c+1:c+3)=f1(size_drum11+1,:);
                temp_ca(:,c+1:c+3)=f1(size_drum11+2,:);
                temp_n(:,c+1:c+3)=f1(size_drum11+3,:);
                temp_o(:,c+1:c+3)=f1(size_drum11+4,:);
                
                count=count+1;
                %%to keep tab of uncommon backbones
                temp_c1(count,:)=f1(size_drum11+1,:);
                temp_ca1(count,:)=f1(size_drum11+2,:);
                temp_n1(count,:)=f1(size_drum11+3,:);
                temp_o1(count,:)=f1(size_drum11+4,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
    end
    
    
    num1=o1;
    num_always=num1;
    clear drum PROB1 ENE1 temp_n1 temp_ca1 temp_c1 temp_o1;
    %pos_prob = PROB1; pos_ener = ENE1;
    drum = temp11; %PROB1=pos_prob; ENE1=pos_ener;
    clear temp11; %pos_prob pos_ener;
    %end
    
end



clear drum11; clear pos_probability; clear pos_energy;
drum11 = drum; pos_probability = pos_prob; pos_energy = pos_ener;
clear drum pos_prob pos_ener temp11 PROB1 ENE1;
clear temp_n temp_ca temp_c temp_o;
temp_c(1,:)=drum11(y1+1,:);
temp_ca(1,:)=drum11(y1+2,:);
temp_n(1,:)=drum11(y1+3,:);
temp_o(1,:)=drum11(y1+4,:);
 
Num_pos_tot(i,3)=num1;
%clear drum;

y1=l; y2=num1;

delete pos1_*.pdb; clear t; clear t1; clear des; clear CC;clear CC1; %delete coor1_*.pdb;
clear iNum iRes rn r1 PR EN TO Begin;
disp(num1)
fname=sprintf('chkpoint.mat');
    save(fname);

end

end


