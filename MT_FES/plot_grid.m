function plot_grid(bin)
%%%%%%%%% Plotting the gridded data for FES %%%%%%%%%%%%%%%%%%%%%%%%
data=load('rncoord_p_pl.txt');

rn_p_pl=[data(:,1) data(:,2) data(:,3) data(:,4)];

xmin=10.2000;
xmax=20.4000;
ymin=14.0000;
ymax=27.8000;

bin=0.555;

rn_p_pl(:,3)=round(rn_p_pl(:,3));
rn_p_pl(:,4)=round(rn_p_pl(:,4));

xx = xmin:bin:xmax;
yy = ymin:bin:ymax;

P=max(rn_p_pl(:,3));
PL=max(rn_p_pl(:,4));

data_P_lowest=zeros((size(xx,2)*size(yy,2)),3);
data_PL_lowest=zeros((size(xx,2)*size(yy,2)),3);

n=0;
for i=1:size(xx,2)
    for j=1:size(yy,2)

        n=n+1;
        data_E = rn_p_pl( abs(rn_p_pl(:,1)-xx(1,i))<=bin/2 & abs(rn_p_pl(:,2)-yy(1,j))<=bin/2,: );

        if isempty(data_E)
             data_P_lowest(n,:)=[xx(1,i) yy(1,j) P];
             data_PL_lowest(n,:)=[xx(1,i) yy(1,j) PL];

        elseif ~isempty(data_E)
            [a,b]=size(data_E);
            clear temp_listP temp_listPL;
            temp_listP = sortrows(data_E,3);
            temp_listPL = sortrows(data_E,4);

           m1=1;

            for j1=1:m1
                data_P_lowest(n+j1-1,:)=[xx(1,i) yy(1,j) temp_listP(j1,3)];
                data_PL_lowest(n+j1-1,:)=[xx(1,i) yy(1,j) temp_listPL(j1,4)];
            end
            n=n+m1-1;
        end

    end
end
temp_dP = sortrows(data_P_lowest,3);
temp_dPL= sortrows(data_PL_lowest,3);

dP=[data_P_lowest(:,1) data_P_lowest(:,2) (data_P_lowest(:,3)-temp_dP(1,3))];
dPL=[data_PL_lowest(:,1) data_PL_lowest(:,2) (data_PL_lowest(:,3)-temp_dPL(1,3))];

clear temp_dP temp_dPL;
temp_dP = sortrows(dP,3);
temp_dPL= sortrows(dPL,3);
[n1,n2]=size(temp_dPL);

for i =1:n1
	if temp_dPL(i,3)>70
		temp_dPL(i,3)=70;
	end
end

xa=temp_dPL(1,1); ya=temp_dPL(1,2); za=temp_dPL(1,3);
di=0.001;
for i=1:4
	a=xa+(i*di); b=ya+(i*di); c=za+(i*di);
	temp_dPL(n1+1,1)=a; temp_dPL(n1+1,2)=b; temp_dPL(n1+1,3)=c;
	n1=n1+1;
end

[n1,n2]=size(temp_dP);
for i =1:n1
        if temp_dP(i,3)>70
                temp_dP(i,3)=70;
        end
end

xa=temp_dP(1,1); ya=temp_dP(1,2); za=temp_dP(1,3);
di=0.001;
for i=1:4
        a=xa+(i*di); b=ya+(i*di); c=za+(i*di);
        temp_dP(n1+1,1)=a; temp_dP(n1+1,2)=b; temp_dP(n1+1,3)=c;
        n1=n1+1;
end


clear dP dPL;
dP=sortrows(temp_dP,1);
dPL=sortrows(temp_dPL,1);


dlmwrite('grid_dP.pmf',dP,'delimiter','\t');
dlmwrite('grid_dPL.pmf',dPL,'delimiter','\t');
end
