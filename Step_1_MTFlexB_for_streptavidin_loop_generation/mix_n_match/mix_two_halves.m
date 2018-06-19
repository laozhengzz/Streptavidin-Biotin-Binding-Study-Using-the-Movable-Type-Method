%%%%%% Mix and Match %%%%%%%%%%%

function mix_two_halves(ru1)
num_pdb=0;
fold=1;
ru1=str2num(ru1);
%for ru1=1:4183
for ru=1:1147
clearvars -except fold num_pdb ru ru1;
half=4;

fname_o=sprintf('/mnt/scratch/bansalnu/MTflexb/parallel_10_17_2016/1mk5_ligandbound/opposite/chkpoint_4_%d.mat',ru1);

load(fname_o);

%%% Reverse the coordinates of the pocket %%%%%%%%%
[b1,b2] =size(drum11);
poc_coor1=zeros(size(drum11));
a=1;
for i=1:half
    %d=nr(i,1);
    bb=strtrim(name_res(i,:));
%     if (strncmp(bb,'GLY',3)==1) 
%         d=4;
%     else d=5;
%     end
    d=nr(i,1);
    clear aa;
    aa = drum11(a:(a+d-1),:);
    a1=aa(1,:);
    a2=aa(3,:);
    aa(1,:)=a2; aa(3,:)=a1;
    clear a1 a2
    poc_coor1((b1-d+1):b1,:)=aa;
    a=a+d;
    b1=b1-d;
end
clear drum11;
Num_2ndhalf=b2/3;

%%%%%%%%%%%%% load the first half of the chkpoint and mix and match the
%%%%%%%%%%%%% combinations %%%%%%
load chk.mat;
poc_a1=char(zeros(size_pocket(1,1),3));
poc_r1=char(zeros(size_pocket(1,1),3));
poc_s1=char(zeros(size_pocket(1,1),6));
poc_chain1=char(zeros(size_pocket(1,1),1));
poc_h1=zeros(size_pocket(1,1),1);
%poc_num1=zeros(length(nr),3);
l=0;
for i = 1:num_res
    bb=strtrim(name_res(i,:));
    [n_atoms, name_atoms, d1_atoms, d2_atoms, d3_atoms, hb] = aa_info(bb);
%     if (strncmp(bb,'GLY',3)==1) 
%         n_atoms=4;
%     else n_atoms=5;
%     end
    for k = 1:n_atoms
        l=l+1;
        poc_a1(l,:)=name_atoms(k,:);
        poc_r1(l,:)=name_res(i,:);
        poc_s1(l,:)=seq_res(i,:);
        poc_chain1(l,:)=chain_res(i,:);
        poc_h1(l,1)=hb(k,1);
        s=strtrim(name_atoms(k,:));
        if (i == half) && ((strncmp(s,'C',1)==1) && (length(s)==1))
            C_pos=l;
        end
    end
end
l11=l;
fname=sprintf('/mnt/scratch/bansalnu/MTflexb/parallel_10_17_2016/1mk5_ligandbound/chkpoint_4_%d.mat',ru);
load(fname);
[b1,b2]=size(drum11);
Num_1sthalf=b2/3;
l=l11;
Final=zeros(l,100000);

f=0;
for i = 1 : Num_2ndhalf
    d=3*(i-1);
    d11=poc_coor1(1,d+1:d+3);  %% N of the new residue
    d12=poc_coor1(2,d+1:d+3);  %% Ca of the new residue
    for j = 1 : Num_1sthalf
        e=3*(j-1);
        d21=drum11(C_pos,e+1:e+3);  %% C of the previous residue
        d22=drum11(C_pos+1,e+1:e+3); %% O of the previous residue
        d23=drum11(C_pos-1,e+1:e+3); %% Ca of the previous residue

        a1=sqrt((d11(1,1)-d21(1,1))^2 + (d11(1,2)-d21(1,2))^2 + (d11(1,3)-d21(1,3))^2); %N-C bond distance
        a2=sqrt((d11(1,1)-d22(1,1))^2 + (d11(1,2)-d22(1,2))^2 + (d11(1,3)-d22(1,3))^2); %N-O angle distance
        a3=sqrt((d11(1,1)-d23(1,1))^2 + (d11(1,2)-d23(1,2))^2 + (d11(1,3)-d23(1,3))^2); %N-Ca angle distance
        a4=sqrt((d12(1,1)-d22(1,1))^2 + (d12(1,2)-d22(1,2))^2 + (d12(1,3)-d22(1,3))^2); %Ca-C angle distance
        a5=sqrt((d12(1,1)-d21(1,1))^2 + (d12(1,2)-d21(1,2))^2 + (d12(1,3)-d21(1,3))^2); %C-Ca angle distance

        %out = check_pocket_combine(poc_h1,poc_a1,poc_r1,poc_s1,poc_chain1,drum11(:,e+1:e+3),poc_coor1(:,d+1:d+3));

        if (a1 < 1.7) && (a2 >= 2.0) && (a3 >= 2.0) && (a4 >= 2.0) && (a5 >= 2.0)
            out = check_pocket_combine(poc_a1, poc_r1, poc_s1,poc_chain1, poc_h1,drum11(:,e+1:e+3),poc_coor1(:,d+1:d+3));
            if (out == 0)
                f=f+1;
                f1=3*(f-1);
                Final(:,f1+1:f1+3)=[drum11(:,e+1:e+3);poc_coor1(:,d+1:d+3)];
            end
        end
    end
end

%if f > 0
%  new=sprintf('%d',fold);
%  mkdir(new);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% saving everything in pdb format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:f
        r=3*(m-1);
        dummy=Final(:,r+1:r+3); %pos_prob=pos_probability(:,m); pos_ener=pos_energy(:,m);prob1=1; ener1=0;
        num_pdb=num_pdb+1
        for zz=1:l
            ss1=poc_a1(zz,:);
            bb11=poc_r1(zz,:);
            ww=poc_s1(zz,:);
            cc=poc_chain1(zz,:);
            yield=sprintf('%1.4s %6d  %2.4s %3.3s %1.1s%1.4s     %7.3f %7.3f %7.3f','ATOM',zz,ss1,bb11,cc,ww,dummy(zz,1),dummy(zz,2),dummy(zz,3));
            %yield=sprintf('%1.4s %6d %2.4s %3.3s %3.2s %3.4s    %6.3f %7.3f %7.3f','ATOM',zz,ss1,bb11,cc,ww,dummy(zz,1),dummy(zz,2),dummy(zz,3));
            mycell={yield};
            t1=sprintf('pos_%d.pdb',num_pdb);
            %/mnt/home/bansalnu/drug_design/jobs_final/total_package/dlmcell(t1,mycell,'-a',' ');
            dlmcell(t1,mycell,'-a',' ');
        end

    end
 %  if f > 0
%	new=sprintf('%d',fold);
%	cd(new);
%	movefile ../pos*.pdb .
%	cd ../
%	fold=fold+1;
 %  end
end
end
