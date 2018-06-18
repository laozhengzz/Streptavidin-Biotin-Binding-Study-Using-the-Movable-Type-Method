clear all;
clc;
bin=0.5;

%%%% Crystal str. X and Y
%Closed
xc=18.6; yc=15.2;

%Open
xo=12.9; yo=23.4;

RcutoffPP=6;
RcutoffPL=6;

data=load('rn_coord_alldata_docked.txt');

E_P = [data(:,1) data(:,2) data(:,3)+data(:,4)./(RcutoffPP*10) + (data(:,5) + data(:,6))];
E_PL =[data(:,1) data(:,2) (data(:,3)+data(:,4)./(RcutoffPP*10) + data(:,9)./(RcutoffPL)  + data(:,10) + data(:,11)*-0.03)];


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

delG_open_binding = Open_PL - Open_P

delG_close_binding = Close_PL - Close_P


