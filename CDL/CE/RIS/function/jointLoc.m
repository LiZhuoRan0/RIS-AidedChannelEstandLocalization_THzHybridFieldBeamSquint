function UE_pos_jointest = jointLoc(angle_BS_est,angle_RIS_est, BS_pos, RIS_pos)
%% obtain the UE's location by the location of the BS and the RIS, the estimated angle from the UE to the BS, and the estimated angle from the UE to the RIS
% inputs:
%       angle_BS_est                =      estimated angle from UE to BS      
%       angle_RIS_est               =      estimated angle from UE to RIS  
%       BS_pos                      =      BS location
%       RIS_pos                     =      RIS location

% outputs:
%     UE_pos_jointest               =      estimated UE's location

%%
k_UE2BS = tan(angle_BS_est);
k_UE2RIS = tan(angle_RIS_est);

k1 = k_UE2BS;
k2 = k_UE2RIS;
b1 = BS_pos.y-k1*BS_pos.x;
b2 = RIS_pos.y-k2*RIS_pos.x;

UE_pos_jointest.x = (b1-b2)/(k2-k1);
UE_pos_jointest.y = (b1*k2-b2*k1)/(k2-k1);
end