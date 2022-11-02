%% Step 3 Part 2 - Sensitivity Analysis

% Sensitivity analysis on PO flushing designs for different
% sediment yield scenarios

%%
addpath(genpath('DCASCADE_3S-Matlab_code'))
 
load('output_Borg_result.mat')

%check all valid sediment yield scenarios
sed_input_param_comb = [50,0.1,1;55,0.1,1;60,0.1,1;40,0.075,1;45,0.075,1;50,0.075,1;55,0.075,1;60,0.075,1;35,0.05,1;40,0.05,1;45,0.05,1;50,0.05,1;55,0.05,1;25,0.025,1;30,0.025,1;35,0.025,1];
  
%prepare output
output_sensitivity = cell(3,4);
output_sensitivity(:,1:3) = output_Borg_result(:,1:3);

%run for each tr.cap equations
for se = 1:3
    
    theta_param_comb = output_Borg_result{se,3};

    eq_output_sensitivity = cell(size(theta_param_comb,1)+1,size(sed_input_param_comb,1) +1);

    %run for each PO flushing configuration
    for st=1:size(theta_param_comb,1)

        theta_param_comb(st,1) = round(theta_param_comb(st,1));        
        theta_param_comb(st,2) = ceil(theta_param_comb(st,2)); 
        theta_param_comb(st,4) = round(theta_param_comb(st,4));
        if theta_param_comb(st,4) >= 3.9 % if frequency is more then 3 years, flushing is not done
            theta_param_comb(st,4) = 100;
        end

        theta_param_comb(st,5) = round(theta_param_comb(st,5));

        eq_output_sensitivity{st+1,1} = theta_param_comb(st,:);
        
        %run for each sediment yield scenario
        for si=1:size(sed_input_param_comb,1) 

            eq_output_sensitivity{1,si+1} = sed_input_param_comb(si,:);

            switch se
                case 1
                    [JJ] = DCASCADE_3S_sens_ATK(sed_input_param_comb(si,:),theta_param_comb(st,:) );
                case 2
                    [JJ] = DCASCADE_3S_sens_EH(sed_input_param_comb(si,:),theta_param_comb(st,:) );
                case 3
                    [JJ] = DCASCADE_3S_sens_EHW(sed_input_param_comb(si,:),theta_param_comb(st,:) );
            end

            eq_output_sensitivity{st+1,si+1} = JJ;

        end

    end

    output_sensitivity(se,4) = eq_output_sensitivity;
    
end

% save final results
save('Step_3_Reservoir_sediment_management/Part_2_sensitivity_analysis/Part_2_Result/output_sensitivity','output_sensitivity', '-v7.3');
