clear all;
clc;
bin = 0.500;
data=load('rncoord_p_pl.txt');

format longG

rn_p_pl=[data(:,1) data(:,2) data(:,3) data(:,4)];

xmin=10.2000;
xmax=20.4000;
ymin=14.0000;
ymax=27.8000;

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
	    
%            if length(temp_listP) >= 20
 %               m1=20;
 %           elseif (length(temp_listP) > 5) && (length(temp_listP) < 20)
  %              m1=5;
  %          else m1=1;
   %         end
                
            for j1=1:m1
                data_P_lowest(n+j1-1,:)=[xx(1,i) yy(1,j) temp_listP(j1,3)];
                data_PL_lowest(n+j1-1,:)=[xx(1,i) yy(1,j) temp_listPL(j1,4)];
            end
            n=n+m1-1;
        end
        
    end
end


%%%%%%% Find the delG_crystal b/w apo and holo structures %%%%%
%%%% Crystal str. X and Y grids. Note that values have shifted a bit from the exact crystal structure value because of the gridding of the data.
%%% Search was done mainly within the grid of the crystal closed and open strcuture. 

E_P=data_P_lowest;
E_PL=data_PL_lowest;

%Closed grid
xc=18.6; yc=15.4;

%Open grid
xo=13.0; yo=23.8;

%%%%%%% Close like %%%%%%%%%%%%
%%%%%% Protein fe
close = E_P((E_P(:,1) >= (xc - bin) & E_P(:,1) <= (xc + bin)) & (E_P(:,2) >= (yc - bin) & E_P(:,2) <= (yc + bin)),:);
val = sortrows(close,3);
Close_P=val(1,3);

%%%%%% Protein-ligand fe
clear close val;
close = E_PL((E_PL(:,1) >= (xc - bin) & E_PL(:,1) <= (xc + bin)) & (E_PL(:,2) >= (yc - bin) & E_PL(:,2) <= (yc + bin)),:);
val = sortrows(close,3);
Close_PL=val(1,3);

%%%%%%% Open like %%%%%%%%%%%%
%%%%%% Protein fe
open = E_P((E_P(:,1) >= (xo - bin) & E_P(:,1) <= (xo + bin)) & (E_P(:,2) >= (yo - bin) & E_P(:,2) <= (yo + bin)),:);
val = sortrows(open,3);
Open_P=val(1,3);

%%%%%% Protein-ligand fe
clear open val;
%%%%%% Protein fe
open = E_PL((E_PL(:,1) >= (xo - bin) & E_PL(:,1) <= (xo + bin)) & (E_PL(:,2) >= (yo - bin) & E_PL(:,2) <= (yo + bin)),:);
val = sortrows(open,3);
Open_PL=val(1,3);

delG_close_to_open_apo = Close_P - Open_P
delG_close_to_open_holo = Close_PL - Open_PL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

