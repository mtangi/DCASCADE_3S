function [ v_sed, Qtr_cap ] = velocity_EH( Slope_reach , Wac_reach , v_reach , h_reach , minvel, phi , psi )
%VELOCITY_EH returns the velocity of the sediment (in m/s) for each sediment
%class for each reach using the Engelund Hansen equations (1967)
%
%OUTPUTS:
% 
% v_sed: [cxn] matrix reporting the velocity for each sediment class c for
%        each reach n [m/s]
%% references
%Engelund, F., and E. Hansen (1967), A Monograph on Sediment Transport in Alluvial Streams, Tekniskforlag, Copenhagen.


%% active layer definition

%active layer as 10% of the water depth 
L_a = 0.1.*h_reach ; %characteristic vertical length scale for transport.

%alternative: active layer as 2*D90 (Parker, 2008)
%L_a = 2*D_finder_3(Fi_r_reach, 90 );

%% sediment velocity with fractional trasport capacity

%  by measuring the trasport capacity for the single sed.classes
%  indipendenty, we obtain different values of sed. velocity

rho_s = 2650; % sediment densit [kg/m^3]
rho_w = 1000; % water density [kg/m^3]
g = 9.81;
dmi = 2.^(-psi)'./1000;

%friction factor
C = (2*g.*Slope_reach.*h_reach)./(v_reach).^2;
%dimensionless shear stress
tauEH = (Slope_reach.*h_reach)./((rho_s/rho_w-1).*dmi);
%dimensionless transport capacity
qEH = 0.05./C.* (tauEH).^(5/2);

%dimensionful transport capacity 
qEH_dim = qEH.*sqrt((rho_s./rho_w-1)*g*(dmi).^3); %m2/s
QS_kg = qEH_dim.*Wac_reach.*rho_s; %kg/s
QS_EH = QS_kg./rho_s; %m3/s

%calculate velocity
v_sed = max( QS_EH./( Wac_reach .* L_a .* (1-phi) ) , minvel);

v_sed(:,L_a==0) = minvel;

end