function [ Qtr_cap, pci ] = flushing_tr_cap(Fi_r_reach , Slope_reach , Q_flushing, n_Man_reach , psi , par_method)

%flushing_tr_cap returns the value of the transport capacity (in Kg/s) 
%for each sediment class in the reach measured using the function for
%flushing transport capacity calculation of IRTCES(1985) 


%% references
% Atkinson, E. "The feasibility of flushing sediment from reservoirs."(1996). page 5 equation 2
       
%% define flushing width

Wac_FL = 12.8.* Q_flushing.^0.5;

%% Transport capacity using bed material fraction approach (BMF, Molinas and Wu, 2000)

dmi_finer = 2.^(-psi)'; %sediment classes diameter (mm)

rho_s = 2650; % sediment densit [kg/m^3]
%rho_w = 1000; % water density [kg/m^3]
%g = 9.81;

PS = (dmi_finer<0.016).*1600 + and(dmi_finer>0.016,dmi_finer<0.1).*650 + (dmi_finer>=0.1).*300;

QS_FL_ts = PS .* Q_flushing.^1.6 .* Slope_reach.^1.2 ./ Wac_FL.^0.6; %t/s
QS_FL = QS_FL_ts.*1000./rho_s; %m3/s

Qtr_cap = QS_FL.*Fi_r_reach;

end
