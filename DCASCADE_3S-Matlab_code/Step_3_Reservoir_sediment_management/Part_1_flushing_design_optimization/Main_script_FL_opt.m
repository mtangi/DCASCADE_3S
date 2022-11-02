%% Step 3 Part 1 - Flushing design optimization

% Coordinated sediment flushing operations are designed 
% using the BORG multi-objective evolutionary algorithm to explore interesting trade-offs 
% between the two conflicting objectives of energy generation and outlet sediment delivery 

% The BORG algorithm for Matlab can be requested at http://borgmoea.org/

%Objectives: 

% J1 = hydropowers (mean GWh/yr )
% J2 = sediment trasport (Mean delivery to the outlet Mt/yr)

% Parameters:

% theta(1) = duration_fl  
% theta(2) = Fl_param.startmonth
% theta(3) .* max(Q(:,[DamDatabase_active(d).reach_UP]),[],1) = minflow
% theta(4) = annfreq
% theta(5) = fllag

% Flushing sediment transport equation
% 1: Atkinson (1996)
% 2: Engelund & Hansen (1967)
% 3: Engelund & Hansen (1967) w/t width correction equation in Atkinson (1996)

%%
addpath(genpath('DCASCADE_3S-Matlab_code'))

nvars = 5 ;
nobjs = 2 ;
NFE   = 20000 ; 
LB    = [1 0 0 1 0]; %lower boundaries for decision variables
UB    = [15 12 1 4 60]; %upper boundaries for decision variables
eps_dom = [ 1 , 0.01 ];

%initialize output cell structure
output_Borg_results = cell(3,3);
output_Borg_results{1,1} = 'flushing_eq';
output_Borg_results{2,1} = 'EH_eq';
output_Borg_results{3,1} = 'EH_eq_changeWac';

%run for each flushing sediment transport equation
for i = 1:3

    switch i
        case 1
            [theta, Jopt,runtime] = borg(nvars, nobjs, 0, @DCASCADE_3S_opt_ATK, NFE, LB, UB, eps_dom);
        case 2
            [theta, Jopt,runtime] = borg(nvars, nobjs, 0, @DCASCADE_3S_opt_EH, NFE, LB, UB, eps_dom);
        case 3
            [theta, Jopt,runtime] = borg(nvars, nobjs, 0, @DCASCADE_3S_opt_EHW, NFE, LB, UB, eps_dom);
    end
    
   output_Borg_results{i,2} = Jopt;
   output_Borg_results{i,3} = theta;
    
end

% save final results
save('Step_3_Reservoir_sediment_management/Part_1_flushing_design_optimization/Part_1_Result/output_Borg_result','output_Borg_result', '-v7.3');
