%% Step 2 - Baseline sediment budget definition

%  D-CASCADE is applied on the 3S basin for three different portfolios of dam siting, 
%  with a special focus on downstream dams (which would trap most of the basin?s sediment). 
%  These portfolios represent alternative dam development configurations for the lower 3S basin.

%%

addpath(genpath('DCASCADE_3S-Matlab_code'))

load('DCASCADE_3S-Matlab_code/Step_2_Reservoir_siting/Step_2_Input/sed_budget_acceptable.mat')

data_output_allcomb = cell(length(sed_budget_acceptable),3);

for s=1:length(sed_budget_acceptable)
    
    [data_output,dam_output] = DCASCADE_3S_damimpact( sed_budget_acceptable(:,s) );
    
    data_output_allcomb{s,1} = data_output;
    data_output_allcomb{s,2} = dam_output;
    data_output_allcomb{s,3} = sed_budget_acceptable(:,s);
    
end

% save final results
save('Step_2_Result/output_damimpact','output_damimpact', '-v7.3');
