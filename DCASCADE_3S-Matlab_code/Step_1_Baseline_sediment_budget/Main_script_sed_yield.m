%% Step 1 - Baseline sediment budget definition

% We apply the basin-scale sediment connectivity model D-CASCADE 
% on the entire 3S river network without reservoirs, to quantify
% reach sediment budgets and transport rates for the unimpeded river system and 
% compensate for the lack of distributed field data. 

% Multiple scenarios with different input catchment sediment yields are analyzed 
% to account for the uncertainty in the model initialization. 

% Using quantitative field observations and literature estimations regarding
% annual sediment transport at the 3S river confluence, we select several scenarios
% to serve as a baseline for the simulations with reservoirs.

%% 

addpath(genpath('DCASCADE_3S-Matlab_code'))

%define sediment yield classes
sed_yield = 15:5:60;
sed_D50 = [1 0.75 0.5 0.25 0.1 0.075 0.05 0.025 0.01 0.0075 0.005];

%find all combinations of the vector
% Use meshgrid to generate all combinations
[A, B] = meshgrid(sed_yield,sed_D50);
A = A(:);
B = B(:);
sed_input_param_comb = [A, B];

%or if you have the Deep Learing Toolbox
%sed_input_param_comb = combvec(sed_yield, sed_D50);

%initialize output elements
data_output_allcomb = cell(length(sed_input_param_comb),2);

% run D-CASCADE for each scenario
for s=1:length(sed_input_param_comb)
    
    [data_output] = DCASCADE_3S_sedyield(sed_input_param_comb(s,:) );
    
    data_output_allcomb{s,1} = data_output;
    data_output_allcomb{s,2} = sed_input_param_comb(s,:);
    
end

% save final results
save('Step_1_Result/output_baseline_sed_budget','data_output_allcomb', '-v7.3');
