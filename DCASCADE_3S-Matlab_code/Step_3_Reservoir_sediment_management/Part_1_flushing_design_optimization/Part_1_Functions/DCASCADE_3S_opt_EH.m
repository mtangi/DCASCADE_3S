function [ JJ] = DCASCADE_3S_opt_EH(theta)
%% CASCADE
%
% INPUT :
%
% AggData        = Struct defining the features of the network reaches
% Network        = 1x1 struct containing for each node info on upstream and downstream nodes
%
%----
% OUTPUT: 
%
% Qbi_tr         = 
%
% QB_tr          = 
%
%CASCADE_3_newvel calculate first the sed.distribution for all reaches and
%then the velocity and the sediment transfer - to do so requires a second reach loop.

%% load data

% Load network data
load('network_data_3S.mat')

Lngt = [ReachData.Length];
n_Man = [ReachData.n];
Wac = [ReachData.Wac];
Slope = [ReachData.Slope];

outlet = find([ReachData.FromN] == [ReachData.ToN]);
Node_el = [[ReachData.el_FN] [ReachData(outlet).el_TN]];

% Load initialization data for CASCADE
load('Q_data_3Sborg.mat')
timescale = length(Q);
Q_original = Q;

% Load DamDatabase
load('DamDatabase_3S.mat')

% load input sediment load
load('Qbi_input_3Sborg_scenarios.mat')


%% round theta

%{
    theta(1) = duration  
    theta(2) = Fl_param.startmonth
    theta(3) .* max(Q(:,[DamDatabase_active(d).reach_UP]),[],1) = minflow
    theta(4) = annfreq
    theta(5) = fllag
%}

theta(1) = round(theta(1)); %duration
theta(2) = ceil(theta(2)); 
theta(4) = round(theta(4));
if theta(4) >= 3.9 % if frequency is more then 3 years, flushing is not done
    theta(4) = 100;
end

theta(5) = round(theta(5));

%% initialize dam variables
load('startstate_3Sborg.mat')

ResVolume = zeros(timescale, length(DamDatabase_active));
ResVolume(1:2,:) = repmat(ResVolume_start,2,1);

release = zeros(timescale, length(DamDatabase_active));
release(1,:) = release_start;

Qbi_tr = Qbi_tr_start;
%Qbi_tr(2,:,:) = Qbi_tr_start;

Qbi_dep = Qbi_dep_start;
%Qbi_dep(2,:,:) = Qbi_dep_start;

Qbi_mob  = zeros(length(ReachData),length(psi_3S));

Q_out = zeros( 1, length(psi_3S));

%% initialize variables
% sediment velocity parameters
phi = 0.4; 
minvel = min([ReachData.Slope]); %0.00001;

%% initialize flushing parameters

%damushing(t,d,z) contains in z=1 the stage of flushing for each
%timestep t for each dam d, and in z=2 the number of days in that stage.

dam_flushing = zeros(timescale, length(DamDatabase_active) , 2);
dam_flushing(1,:,1) = ones(1,length(DamDatabase_active));
dam_flushing(1,:,2) = ones(1,length(DamDatabase_active));

%% set damushing parameters

for d=1:length(DamDatabase_active)
    
    ORparameters_active(d).Fl_param.duration = theta(1);
    
    ORparameters_active(d).Fl_param.startmonth = theta(2);
    ORparameters_active(d).Fl_param.minflow = theta(3).* max(Q(:,[DamDatabase_active(d).reach_UP]),[],1) ;
        
    ORparameters_active(d).Fl_param.annfreq = theta(4);

    if d==2
       ORparameters_active(d).Fl_param.fllag = theta(5); 
    else
       ORparameters_active(d).Fl_param.fllag = 0; 
    end
    
end

clear mnth
 flushing_yes = [nan, 0, nan];
% 
%% Routing scheme
% initialize new variables 
% Routing scheme  
for t = 2:timescale-1
      
    %variables initialization
    %Qbi_tr(2,:,:) = zeros( size(Qbi_tr(1,:,:)) );
    %Qbi_dep(2,:,:) = zeros( size(Qbi_dep(1,:,:)) );
    Qbi_mob = zeros( size(Qbi_mob) );
    
    % change reservoir FSL according to reservoir sediment storage
    [FSL_ResVolume] = changeFSLResVolume(DamDatabase_active,  Qbi_tr ,  Qbi_dep , flooded_reaches_full,phi);     

    %calculate dam mass balance and dam release with flushing, for the LSS2_II dam which must syncronize with upstream reservoirs US1
    if ~isnan(flushing_yes(2))
        
        if t>ORparameters_active(2).Fl_param.fllag && flushing_yes(2) == 0 && dam_flushing(t-1,2,1) >=0 && dam_flushing(t-ORparameters_active(2).Fl_param.fllag,1,1)<=-1  
            flushing_yes(2)=1;
        else
            flushing_yes(2)=0;
        end
    
    end

    %calculate dam mass balance and dam release with flushing
    [release(t,:), ResVolume(t+1,:), Q(t,:) , dam_flushing(t,:,:)] = dam_mass_balance_flushing(DamDatabase_active, ORparameters_active, Network,   ResVolume(t,:), Q_original(t,:), Q_original(t-1,[DamDatabase_active.reach_UP]), dates_Q(:,t),  FSL_ResVolume,  squeeze(dam_flushing(t-1,:,:)) , WL_target(t,:), flushing_yes);
        
    %calculate reach features for the flooded reaches
    [ResInundNodes , Wreach, Energy_slopereach,  hreach, avgow_vreach ] = dam_features_correction( DamDatabase_active , ReachData, Network, ResVolume(t,:) , Slope , Node_el);
    
    %extract values of ReservoirHydraulics
    f_reach = cell2mat(ResInundNodes');    
    
    EnergySlope = Slope; %in all reaches but the flooded ones, energy slope is equal to the slope 
    EnergySlope(f_reach) = min( cell2mat(Energy_slopereach'), Slope(f_reach) );
    
    Wac(f_reach) = max( cell2mat(Wreach'), [ReachData(f_reach).Wac] );
    
    %calculate new water dept for all reaches
    %Manning
    h = (Q(t,:).*n_Man./(Wac.*sqrt( EnergySlope ))).^(3/5);
    h(f_reach) = max( cell2mat(hreach'), h(f_reach));

    v = 1./n_Man.*h.^(2/3).*sqrt( EnergySlope ); 
    v(f_reach) = min( cell2mat(avgow_vreach') , v(f_reach));
        
    for n = NH_dams

        % extract the deposit layer of the reach from the relative cell in the previous timestep
        V_dep_old = Qbi_dep(n,:);

        %%% 1) extract the deposit layer from the storage matrix and load the incoming cascades

        Qbi_incoming = Qbi_tr(n,:) +  Qbi_inputFL{t}(n,:)  ; 

        %%% 2) find cascades to be included into the active layer according to the limit V_lim_tot, and use the cumulative GSD to compute tr_cap

        % find total sediment volume GSD
        Fi_r_reach =  ( V_dep_old + Qbi_incoming ) ./  sum( ( V_dep_old + Qbi_incoming ) ,'all') ; %i find the GSD of the active layer, for the transport capacity calculation
        Fi_r_reach(isnan(Fi_r_reach)) = 0;  %if V_act is empty, i put Fi_r equal to 0 for all classes
        
        %calculate transport capacity using the Fi of the active layer, the resulting tr_cap is in m3/s and is converted in m3/day
        % in case the active layer is empty, the tr_cap is equal to 0
        if sum(Fi_r_reach)==0
            tr_cap = zeros(size(psi_3S));
        else
            tr_cap = Engelund_Hansen_tr_cap(Fi_r_reach' , EnergySlope(n) , Wac(n), v(n) , h(n) ,psi_3S)' .* 24.*60.*60;
        end 

        %%% 3) Deposit the cascades in the active layer until the volume mobilized for each class is equal to the tr_cap
        V_mob = min( V_dep_old + Qbi_incoming , tr_cap);

        V_dep = max(0, V_dep_old + Qbi_incoming - V_mob);

        % save the deposit volume 
        Qbi_dep(n,:) = V_dep;

        % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream
        Qbi_mob(n,:) = V_mob;

    %(end of the reach loop)
    end
    %%% 5) Move the mobilized volumes to the destination reaches according to the sediment velocity

    %loop for all reaches, now that i have the Fi_r and thus can compute transfer rates for all reaches
    clear Qbi_tr_t

    [ v_sed ] = velocity_EH( EnergySlope , Wac , v , h , minvel, phi ,psi_3S) ;

    Qbi_tr = zeros( size(Qbi_tr) );
   
    for n = NH_dams

        V_mob = Qbi_mob(n,:);

        if sum(V_mob,'all') >0

            [Qbi_tr_t, Q_out_t  ] = sed_transfer( V_mob , n , v_sed.*(60*60*24) , Lngt, Network );

            Qbi_tr = Qbi_tr + Qbi_tr_t;
            Q_out = Q_out + Q_out_t;
        end
    end    
        
    %add the external sediment contribution arriving directly to the outlet
    Q_out = Q_out + Q_out_extraFL(t,:);
        
%end of the time loop  
end
%% calculate tot sediment delivery

tot_sed_del = sum(Q_out,'all')*2.6./1000000 /(timescale/365);

%% calculate energy production

dam_heigth = zeros(timescale,size(release,2));
energy_release = zeros(size(release));

for d=1:length(DamDatabase_active)
    dam_heigth(:,d) = Reservoir_LSWConversion( [ResVolume(1:timescale,d)] ,DamDatabase_active(d).WL_table, 3, 1);
    energy_release(:,d) = min(release(:,d) , DamDatabase_active(d).design_discharge) ;
    energy_release( or(dam_flushing(:,d,1) == -1, dam_flushing(:,d,1) == -2),d) = 0;   
end

power_output = 0.9 .* 1000 .* min(energy_release, [DamDatabase_active.design_discharge]) .* 9.81 .* dam_heigth ./10^9; %power generation GW

tot_energy_output = sum(power_output.*24,'all')/(timescale/365); %energy generation - GWh/y

%% save objective results

JJ = [ -tot_energy_output , -tot_sed_del];

%% save results in txt file

% a = [datenum(datetime('now')) , (theta) , JJ ];
% dlmwrite('Results/DCASCADE_borgrun/cascaderuns.txt',a,'-append','coffset',1);

end