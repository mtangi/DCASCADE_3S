function [release_t_dam, minrelease_new, dam_flushing_post_dam] = flushing_relase(DamDatabase_dam, ORparameters_dam , dam_flushing_pre_dam , date_timestep, ResVolume_pre_dam , Q_input_dam , Q_input_estimates_dam, hp_maxrelease , minrelease, flushing_yes_dam)

%FLUSHING_RELEASE determines the release of the dam during flushing
%operation
%Flushing duration is defined a priori.
%Flushing start can be decided no matter the conditions via the input parameter flushing_yes

%% import flushing parameters

month_start = ORparameters_dam.Fl_param.startmonth;
window = ORparameters_dam.Fl_param.flwindow;
months = mod(month_start:month_start+window-1,12);
months(months==0) = 12;

duration_dd = ORparameters_dam.Fl_param.duration_dd;
duration_fl = ORparameters_dam.Fl_param.duration_fl;
duration_ps = ORparameters_dam.Fl_param.duration_ps;

minflow = ORparameters_dam.Fl_param.minflow;
annfreq = ORparameters_dam.Fl_param.annfreq;
deltah_max = ORparameters_dam.Fl_param.deltah_max;

stage_past = dam_flushing_pre_dam(1);
n_days_past = dam_flushing_pre_dam(2);

dam_flushing_post_dam = dam_flushing_pre_dam;

% dam_flushing_stage:
% 0 : no flushing needed
% n : n years without flushing
% -1 : drawdown
% -2 : flushing
% -3 : pause


%% identify stage for the timestep

stage = stage_past;

if stage_past == -3 && n_days_past >= duration_ps %if the number of days from the last failed attempt at flushing exceeds the threshold (pause)...
    stage = annfreq;  % ... try flushing again.
end
        
if stage_past >= 0 %if flushing still needs to be done...
    %...check if the extremes for flushing are respected 
   
    if isnan(flushing_yes_dam) %if flushing_yes=NaN, the flushing is performed authomatically by the function when the conditions are rigth

        check_freq = stage_past >= annfreq; %check if the number of years since last flushing is higher than the annual frequency of flushing;
        check_month = any(date_timestep(2) == months); %check if the timestep falls into the time window for flushing
        check_Q =  Q_input_estimates_dam > minflow; %check if the estimated input diascharge is more than the minimum for flushing
        check_Q_max =  Q_input_estimates_dam < DamDatabase_dam.design_discharge; %check if the estimated input discharge is less than the design discharge

        if any([ check_freq, check_month, check_Q, check_Q_max] == 0) % if any of the check values is not true....      
            % ...flushing is not performed, operations continue as normal

            %if I enter a new year, increase the count of year without flushing
            if and(date_timestep(2) == 1, date_timestep(3) == 1)
                stage = stage_past+1;
            else
                stage = stage_past;
            end

        else % if all checks are passed, start drawdown
            stage = -1;
        end
        
    else %else, the user have specified if flushing must be done
        
        if flushing_yes_dam == 0 % if flushing is not pushed by the user
            if and(date_timestep(2) == 1, date_timestep(3) == 1)
                stage = stage_past+1;
            else
                stage = stage_past;
            end
        else % if all checks are passed, start drawdown
            stage = -1;
        end
    end
end

% if the reservoir is doing flushing operations...

if stage_past == -1 && round(ResVolume_pre_dam,0) == 0 % if drawdown is finished and reservoir is empty...
    stage = -2;  % ...start flushing
elseif stage_past == -1 && n_days_past >= duration_dd %or(n_days_past >= duration_dd, Q_input_dam >= hp_maxrelease + DamDatabase_dam.bottom_discharge) %if the drawdown is taking too long, or the input discharge surpasses the design+bottom discharge ... 
    stage = -3; %...Stop it and pause for a fixed period of time
elseif stage_past == -1 && ResVolume_pre_dam > 0  % if drawdown is continuing and reservoir is still not empty...
    stage = stage_past; % ...continue drawdown
elseif stage_past == -2 && n_days_past >= duration_fl %if the max number of days for flushing is exceeded ...
    stage = 0; %...Start refill
elseif stage_past == -2 && Q_input_dam >= DamDatabase_dam.bottom_discharge %if the input discharge during flushing exceeds the maximum discarge from the bottom outlet...
    stage = -3; %...Stop it and pause for a fixed period of time
elseif stage_past == -3 && ~any(date_timestep(2) == months) %if the dam was waiting after a failed flushing but the time window ended...
    stage = annfreq;  % ... return to normal schedule
end

% release_t = -1 indicates that the normal release schedule is unchanged and flushign is not happening
release_t_dam = -1;
minrelease_new = minrelease;
       
%% define release if flushing is active

switch stage
    
    case -1 %drawdown phase
        
        h_past = Reservoir_LSWConversion( min(ResVolume_pre_dam , DamDatabase_dam.FSL_ResVolume) , DamDatabase_dam.WL_table, 3, 1);
        max_volume_release =  Reservoir_LSWConversion(  max( h_past - deltah_max,0) ,DamDatabase_dam.WL_table, 1, 3);
        maxdrawdown = (ResVolume_pre_dam + Q_input_dam *(24*60*60)/1E6 - max_volume_release)/60/60/24*10^6; %maximum drawdown discharge (Q-input - release) to avoid decrease in elevation = deltah_max
        
        maxrelease =  min(hp_maxrelease + DamDatabase_dam.bottom_discharge);
        opt_release = maxrelease;  %optimal release for the drawdown is the max possible release (limted by the release 
        
        %actual release given the lower of the upper release limits (to avoid negative storage and to avoid excessive level reduction)
        %i also must release at least all incoming discharge  
        release_t_dam = max(Q_input_dam, min( min(maxrelease,maxdrawdown) , opt_release ));
        
    case -2 %flushing phase
        
        release_t_dam = Q_input_dam; %release what enters in the reservoir
        
%     case -3 %refill phase
%         %(the refill phase is conducted as a normal operation phase)
%         minrelease_new = min(Q_input_dam,Fl_mirelease_refill); %minimum discharge during refill (it can be different than the general min discharge to respect seasonality on the hydrology)
end

%% change number of days in the stage

if stage == stage_past % if stage did not change
    n_days_new = n_days_past+1; %continue day count
else %if stage changed
    n_days_new = 1; %start count over
end
    
%% define outputs

dam_flushing_post_dam(1) = stage;
dam_flushing_post_dam(2) = n_days_new;

end

