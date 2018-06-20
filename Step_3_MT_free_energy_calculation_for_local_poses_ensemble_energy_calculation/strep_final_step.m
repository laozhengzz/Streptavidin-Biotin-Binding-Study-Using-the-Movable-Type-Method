clear all
load '/Users/Nupur/Desktop/close_like.txt'
load '/Users/Nupur/Desktop/open_like.txt'

RcutoffPP = 6;
RcutoffPL = 6;
%P_Tor	P_Non	ZE_P_san	ZE_L_san	ZE_P_water	ZE_L_water	COM 	ZE_PL_san	ZE_PL_water

E_P_CL_ind(:,1) = close_like(:,1)+close_like(:,2)./(RcutoffPP*10) + (close_like(:,3) + close_like(:,4));
E_PL_CL_ind(:,1) = close_like(:,1)+close_like(:,2)./(RcutoffPP*10) + close_like(:,7)./(RcutoffPL)  + close_like(:,8) + close_like(:,9)*-0.05;

E_P_OP_ind(:,1) = open_like(:,1)+open_like(:,2)./(RcutoffPP*10) + (open_like(:,3) + open_like(:,4));
E_PL_OP_ind(:,1) = open_like(:,1)+open_like(:,2)./(RcutoffPP*10) + open_like(:,7)./(RcutoffPL)  + open_like(:,8) + open_like(:,9)*-0.05;

temp_list1 = sortrows(E_P_CL_ind(:,1),[1]);
E_P_CL_ind_low(:,1) = temp_list1(1:100,:)-mean(temp_list1(1:100,:));

temp_list2 = sortrows(E_PL_CL_ind(:,1),[1]);
E_PL_CL_ind_low(:,1) = temp_list2(1:100,:)-mean(temp_list2(1:100,:));

temp_list3 = sortrows(E_P_OP_ind(:,1),[1]);
E_P_OP_ind_low(:,1) = temp_list3(1:100,:)-mean(temp_list3(1:100,:));

temp_list4 = sortrows(E_PL_OP_ind(:,1),[1]);
E_PL_OP_ind_low(:,1) = temp_list4(1:100,:)-mean(temp_list4(1:100,:));

Ensemble_E(1,1)=sum(temp_list1(1:100,:).*exp(E_P_CL_ind_low(:,1)./-0.5918)./sum(exp(E_P_CL_ind_low(:,1)./-0.5918)));%+mean(temp_list1(1:100,:));
Ensemble_E(2,1)=sum(temp_list2(1:100,:).*exp(E_PL_CL_ind_low(:,1)./-0.5918)./sum(exp(E_PL_CL_ind_low(:,1)./-0.5918)));%+mean(temp_list2(1:100,:));
Ensemble_E(3,1)=sum(temp_list3(1:100,:).*exp(E_P_OP_ind_low(:,1)./-0.5918)./sum(exp(E_P_OP_ind_low(:,1)./-0.5918)));%+mean(temp_list3(1:100,:));
Ensemble_E(4,1)=sum(temp_list4(1:100,:).*exp(E_PL_OP_ind_low(:,1)./-0.5918)./sum(exp(E_PL_OP_ind_low(:,1)./-0.5918)));%+mean(temp_list4(1:100,:));

State_E_min(1,1) = min( E_P_CL_ind(:,1));
State_E_min(2,1) = min( E_PL_CL_ind(:,1));

State_E_min(3,1) = min( E_P_OP_ind(:,1));
State_E_min(4,1) = min( E_PL_OP_ind(:,1));

[sem1,r1] = find(E_P_CL_ind(:,1) == min( E_P_CL_ind(:,1)));
[sem2,r2] = find(E_PL_CL_ind(:,1) == min( E_PL_CL_ind(:,1)));
[sem3,r3] = find(E_P_OP_ind(:,1) == min( E_P_OP_ind(:,1)));
[sem4,r4] = find(E_PL_OP_ind(:,1) == min( E_PL_OP_ind(:,1)));

Path(1,1) = State_E_min(1,1) - State_E_min(3,1);
Path(2,1) = State_E_min(2,1) - State_E_min(1,1);
Path(3,1) = State_E_min(4,1) - State_E_min(3,1);
Path(4,1) = State_E_min(2,1) - State_E_min(4,1);

Path_Ensemble(1,1) = Ensemble_E(1,1) - Ensemble_E(3,1);
Path_Ensemble(2,1) = Ensemble_E(2,1) - Ensemble_E(1,1);
Path_Ensemble(3,1) = Ensemble_E(4,1) - Ensemble_E(3,1);
Path_Ensemble(4,1) = Ensemble_E(2,1) - Ensemble_E(4,1);










