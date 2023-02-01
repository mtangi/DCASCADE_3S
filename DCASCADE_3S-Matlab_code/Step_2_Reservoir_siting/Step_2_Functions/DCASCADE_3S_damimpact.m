function [data_output,dam_output] = DCASCADE_3S_damimpact(sed_input_param )
%% D-CASCADE
%
%DCASCADE_3S_damimpact run an efficient version of D-CASCADE that do
%not trace deposit layering, to save computational time.
%
% Sediment provenance is still traced

%Dam are included in this version

%% Variables extraction from ReachData 
% Load network data
load('network_data_3S.mat')

% Load Q data
load('Q_data_3S.mat')

%% time initialization 

timescale = length(Q);

%% load old variables
NH = Network.NH;

Wac = [ReachData.Wac];
Slope = [ReachData.Slope];
Lngt = [ReachData.Length];
n_Man = [ReachData.n];

outlet = find([ReachData.FromN] == [ReachData.ToN]);

clear Node_el ; Node_el = [[ReachData.el_FN] [ReachData(outlet).el_TN]];

%% define input sediment load

%define magnitude of the input sed yield , in t/km2/y
tot_yield = sed_input_param(1).*1e6 ; %total yield for the 3S, in t/y
yield_km_y = tot_yield/sum([ReachData.sed_yield]/max([ReachData.sed_yield]).*[ReachData.directAd]); %yield in t/km2/y, for the TVP geomorphic region
sy = num2cell([ReachData.sed_yield]/max([ReachData.sed_yield]) .* yield_km_y);
[ReachData.sed_yield] = sy{:}; %attribute value of flow according to the percentile previously calculated

%define GSD of the input sed yield
Fi_yield = repmat(Fi_extraction( sed_input_param(2)./1000 , 0.8,psi_3S ),[length(NH),1]) ;

%calculatre input sed. yield
[Qbi_input] = initialize_sed_input(ReachData, Fi_yield, Q, dates_Q,psi_3S, 10);

%% dam scenario selection

%WL_target report for each timestep, for each dam the target reseroir
%heigth for the day.

% dam scenario
scenario = sed_input_param(3); %select dam scenario

switch scenario
    case 0; WL_target = zeros(timescale,1);%scenario no dams
    case 1; load('data_LSS2only.mat'); %scenario only LSS2 
    case 2; load('data_fulldam'); %scenario all large dams 
    case 3; load('data_fulldamalt');%scenario all dams - alternative
end

if scenario==0   %scenario no dams
    load('DamDatabase_reduced_full.mat'); 
    p = num2cell([0 0 0 0 0 0 0 ]); 
    [DamDatabase.portfolio] = p{:}; 
    DamDatabase_active = DamDatabase([DamDatabase.portfolio] == 1) ;
    ORparameters_active = ORparameters([DamDatabase.portfolio] == 1) ;
    [ResInundNodes_full ] = dam_features_correction( DamDatabase_active , ReachData, Network, [DamDatabase_active.FSL_ResVolume], Slope(1,:) , Node_el(1,:));
    flooded_reaches_full = ResInundNodes_full;
    %attribute to each flooded reach the ID of the dam it belongs to
    reaches_dam = zeros(size(ReachData))';
    for d = 1:length(DamDatabase_active)   
        reaches_dam(flooded_reaches_full{d} ) = d ;
    end
end

for d = 1:length(DamDatabase_active) 
    DamDatabase_active(d).WL_start = WL_target(1,d); % find target water level
end

%% initialize variables
% sediment velocity parameters
phi = 0.4; 
minvel = min([ReachData.Slope]); %0.00001;

%% initialize dam variables
clear start_vol

start_vol = zeros(1,length(DamDatabase_active));
for d=1:length(DamDatabase_active) 
    if isempty(DamDatabase_active)
        start_vol(d)=0;
    else
        start_vol(d) = Reservoir_LSWConversion(DamDatabase_active(d).WL_start ,DamDatabase_active(d).WL_table, 1, 3);
    end
end

ResVolume = zeros(timescale+1, length(DamDatabase_active));
ResVolume(1:2,:) = repmat(start_vol,2,1);

release = zeros(timescale+1, length(DamDatabase_active));

FSL_ResVolume = zeros(timescale+1, length(DamDatabase_active));
FSL_ResVolume(1:2,:) = repmat([DamDatabase_active.FSL_ResVolume],2,1);

sed_storage = zeros(timescale+1, length(DamDatabase_active));

%operating rule ID
id_OR = 3;

%% Routing scheme
% initialize new variables 

clear Qbi_tr; Qbi_tr  = cell(2, 1);
clear Qbi_mob; Qbi_mob  = cell(1, 1);
clear Qbi_dep; Qbi_dep  = cell(2, 1);
clear Q_out; Q_out = cell(timescale,1);
clear EnergySlope; EnergySlope = zeros(size(Slope));

Qbi_tr{1} = zeros([size(Network.II),length(psi_3S)]);  %Qbi_tr report the sediment mobilized present in the reach AFTER transfer
Qbi_mob{1} = zeros([size(Network.II),length(psi_3S)]); %Qbi_mob report the sediment mobilized present in the reach BEFORE transfer
Qbi_dep{1} = zeros([size(Network.II),length(psi_3S)]);
Qbi_tr{2} = Qbi_tr{1};
Qbi_dep{2} = Qbi_dep{1};
Q_original = Q;

Q_out{1} = zeros( length(Network.II), length(psi_3S));

% Routing scheme  

for t = 2:timescale-1
    
    %variables initialization
    Qbi_tr{2} = zeros( size(Qbi_tr{1}) );
    Qbi_dep{2} = zeros( size(Qbi_dep{1}) );
    Qbi_mob{1} = zeros( size(Qbi_mob{1}) );
    Q_out{t} = zeros( size(Q_out{1}) );
        
    % change reservoir FSL according to reservoir sediment storage
    [FSL_ResVolume(t,:), sed_storage(t,:)] = changeFSLResVolume(DamDatabase_active,  squeeze(sum(Qbi_tr{1},1)) ,  squeeze(sum(Qbi_dep{1},1)) , flooded_reaches_full,phi);     

    %calculate dam mass balance and dam release with flushing
    if scenario ~= 0
       [release(t,:), ResVolume(t+1,:), Q(t,:)] = dam_mass_balance_fast(DamDatabase_active, Network,   ResVolume(t,:), Q_original(t,:), Q_original(t-1,[DamDatabase_active.reach_UP]),  FSL_ResVolume(t,:) , WL_target(t,:));
    end
    
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

    for n = NH

        % extract the deposit layer of the reach from the relative cell in the previous timestep
        V_dep_old = squeeze(Qbi_dep{1}(:,n,:));

        %%% 1) extract the deposit layer from the storage matrix and load the incoming cascades

        Qbi_incoming = squeeze(Qbi_tr{1}(:,n,:)) ;
        Qbi_incoming(n,:) = Qbi_input{t}(n,:);

        %%% 2) find cascades to be included into the active layer according to the limit V_lim_tot, and use the cumulative GSD to compute tr_cap

        % find total sediment volume GSD
        Fi_r_reach =  sum(( V_dep_old + Qbi_incoming ),1) ./  sum( ( V_dep_old + Qbi_incoming ) ,'all') ; %i find the GSD of the active layer, for the transport capacity calculation
        Fi_r_reach(isnan(Fi_r_reach)) = 0;  %if V_act is empty, i put Fi_r equal to 0 for all classes
                
        %calculate transport capacity using the Fi of the active layer, the resulting tr_cap is in m3/s and is converted in m3/day
        % in case the active layer is empty, the tr_cap is equal to 0
        tr_cap = Engelund_Hansen_tr_cap(Fi_r_reach' , EnergySlope(n) , Wac(n), v(n) , h(n) ,psi_3S)' .* 24.*60.*60;

        %%% 3) Deposit the cascades in the active layer until the volume mobilized for each class is equal to the tr_cap
        
        prc = min( tr_cap  ./ sum(( V_dep_old + Qbi_incoming ),1) ,1);
        
        V_mob = [V_dep_old + Qbi_incoming].*prc ;

        V_dep = max(0, V_dep_old + Qbi_incoming - V_mob);

        % save the deposit volume 
        Qbi_dep{2}(:,n,:) = V_dep;
        
        % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream
        Qbi_mob{1}(:,n,:) = V_mob;

    %(end of the reach loop)
    end
    
    %%% 5) Move the mobilized volumes to the destination reaches according to the sediment velocity

    %loop for all reaches, now that i have the Fi_r and thus can compute transfer rates for all reaches
    clear Qbi_tr_t

    [ v_sed ] = velocity_EH( EnergySlope , Wac , v , h , minvel, phi) ;

    for n = NH

        V_mob = squeeze(Qbi_mob{1}(:,n,:));

        if sum(V_mob,'all') >0

            [reach_dest, setout] = sed_transfer_fast( n , v_sed.*(60*60*24) , Lngt, Network );
            
            % i sum the volumes transported in reach n with all the other
            % volumes mobilized by all the other reaches in time t

            % Sum the volumes transported from reach n with all the other
            % volumes mobilized by all the other reaches at time t

            Qbi_tr{t+1} = Qbi_tr{t+1} + single(Qbi_tr_t);
            Q_out{t} =  Q_out{t} + Q_out_t;


            for c=1:length(psi_3S)
                if setout(c) == 0 %if the volume does not leave throught the outlet
                    Qbi_tr{2}(:,reach_dest(c),c)  = Qbi_tr{2}(:,reach_dest(c),c) + V_mob(:,c);
                else 
                    Q_out{t}(:,c) = Q_out{t}(:,c) + V_mob(:,c);
                end
            end

        end
    end    

    
    Qbi_dep{1} = Qbi_dep{2};
    Qbi_tr{1} = Qbi_tr{2};

    
    %end of the timestep skip     
%end of the time loop  
end

%% output processing

outcum_tot = cell2mat(cellfun(@(x) sum(x,'all'), Q_out(1:timescale-1), 'UniformOutput',0));
clear a tot_sed_year D50_year Fi_year
for i=1:timescale/365

    a(i) = sum(outcum_tot(i+(365)*(i-1):365*i-1));
    tot_sed_year(i) = a(i)*2.6./1000000;
    D50_year(i) = D_finder(  sum(cell2mat(cellfun(@(x) sum(x,1), Q_out(i+(365)*(i-1):365*i-1), 'UniformOutput',0)),1)'./a(i) ,psi_3S, 50 )'*1000;
    Fi_year(:,i) = sum(cell2mat(cellfun(@(x) sum(x,1), Q_out(i+(365)*(i-1):365*i-1), 'UniformOutput',0)),1)'./a(i);
end
 
%% calculate  outcum_provenance_river

data_plot_out = Q_out;

river_reach{1,1} = 'Se Kong';
river_reach{1,2} = [find(Network.Upstream.Distance{74} ~= Inf),75];

river_reach{2,1} = 'Se San';
river_reach{2,2} = [find(Network.Upstream.Distance{162} ~= Inf)];

river_reach{3,1} = 'Sre Pok';
river_reach{3,2} = [find(Network.Upstream.Distance{243} ~= Inf)];

river_reach{4,1} = 'Se San - Sre Pok confluence';
river_reach{4,2} = [163:167];

%outcum_provenance contains the mean annual sediment delivery to the outlet
%for each reach, for each sediment class
outcum_provenance = zeros(size(data_plot_out{1}));
data_plot_out{end} = data_plot_out{1};
% 
for n=1:size(data_plot_out{1},1)
    
    outcum_provenance_n_year = zeros(timescale/365,size(data_plot_out{1},2));
    for y=1:timescale/365
        for t=1:365
            outcum_provenance_n_year(y,:) = outcum_provenance_n_year(y,:) + data_plot_out{365*(y-1)+t}(n,:);
        end
    end
    
    outcum_provenance(n,:) = round(mean(outcum_provenance_n_year,1));
end
    
% for t=365:length(Q_out)-1
%     outcum_provenance = outcum_provenance + data_plot_out{t}; 
% end 

%outcum_provenance contains the mean annual sediment delivery to the outlet
%for each reach, for each sediment class
outcum_provenance_river = zeros(size(river_reach,1)+1,size(data_plot_out{1},2) );
for i=1:size(river_reach,1)
    outcum_provenance_river(i,:) = sum(outcum_provenance(river_reach{i,2},:),1);
end
outcum_provenance_river(end,:) = sum(outcum_provenance_river(1:size(river_reach,1),:),1);

%% output struct definition
% data_plot contais the most important D_CASCADE outputs

data_output = cell(1,2);

data_output{1,1} = 'outcum_provenance_river';
data_output{2,1} = 'tot_sed_year';
data_output{3,1} = 'D50_year';
data_output{4,1} = 'Fi_year';
data_output{5,1} = 'Fi_input';

data_output{1,2} = outcum_provenance_river; 
data_output{2,2} = tot_sed_year;
data_output{3,2} = D50_year;
data_output{4,2} = Fi_year;
data_output{5,2} = Fi_extraction( sed_input_param(2)./1000 , 0.8,psi_3S );

if scenario~=0
    dam_output = cell(1,2);

    dam_output{1,2} = DamDatabase_active; 
    dam_output{2,2} = ORparameters_active; 
    dam_output{3,2} = ResVolume;
    dam_output{4,2} = release;
    dam_output{5,2} = FSL_ResVolume;
    dam_output{6,2} = sed_storage; 

    dam_output{1,1} = 'DamDatabase_active'; 
    dam_output{2,1} = 'ORparameters_active';
    dam_output{3,1} = 'Reservoir Volume';
    dam_output{4,1} = 'Release';
    dam_output{5,1} = 'Volume at FSL';
    dam_output{6,1} = 'Reservoir sediment storage';

end

%% all other outputs are included in the extended_output cell variable

% extended_output = cell(1,2);
% 
% extended_output{1,2} = Qbi_tr; 
% extended_output{2,2} = Qbi_mob;
% extended_output{3,2} = Q_out;
% extended_output{4,2} = Qbi_dep;
% extended_output{5,2} = D50_sed;
% extended_output{6,2} = Node_el;
% extended_output{7,2} = Slope;
% 
% extended_output{1,1} = 'Qbi_tr'; 
% extended_output{2,1} = 'Qbi_mob';
% extended_output{3,1} = 'Q_out';
% extended_output{4,1} = 'Qbi_dep';
% extended_output{5,1} = 'D50_AL';
% extended_output{6,1} = 'Node_el';
% extended_output{7,1} = 'Slope';

end