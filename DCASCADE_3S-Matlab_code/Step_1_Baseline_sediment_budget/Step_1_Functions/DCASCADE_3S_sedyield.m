function [data_output] = DCASCADE_3S_sedyield(sed_input_param )
%% DCASCADE_3S_sedyield
%
%DCASCADE_3S_sedyield run an efficient version of D-CASCADE that do
%not trace deposit layering, to save computational time. Sediment provenance is still traced

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

clear Node_el ; Node_el = repmat ([[ReachData.el_FN] [ReachData(outlet).el_TN]], [timescale,1]);

%% define input sediment load

%define magnitude of the input sed yield , in t/km2/y
tot_yield = sed_input_param(1).*1e6 ; %total yield for the 3S, in t/y
yield_km_y = tot_yield/sum([ReachData.sed_yield]/max([ReachData.sed_yield]).*[ReachData.directAd]); %yield in t/km2/y, for the TVP geomorphic region
sy = num2cell([ReachData.sed_yield]/max([ReachData.sed_yield]) .* yield_km_y);
[ReachData.sed_yield] = sy{:}; %attribute value of flow according to the percentile previously calculated

%define GSD of the input sed yield
Fi_yield = repmat(Fi_extraction( sed_input_param(2)./1000 , 0.8,psi_3S ),[length(NH),1]) ;

%define type of split option
split = {'prop'} ;

%calculatre input sed. yield
[Qbi_input] = initialize_sed_input(ReachData, Fi_yield, Q, dates_Q,psi_3S, 10);

%% initialize variables
% sediment velocity parameters
phi = 0.4; 
minvel = min([ReachData.Slope]); %0.00001;

%% Routing scheme
% initialize new variables 

clear Qbi_tr; Qbi_tr  = cell(2, 1);
clear Qbi_mob; Qbi_mob  = cell(1, 1);
clear Qbi_dep; Qbi_dep  = cell(2, 1);
clear Q_out; Q_out = cell(timescale,1);

Qbi_tr{1} = zeros([size(Network.II),length(psi_3S)]);  %Qbi_tr report the sediment mobilized present in the reach AFTER transfer
Qbi_mob{1} = zeros([size(Network.II),length(psi_3S)]); %Qbi_mob report the sediment mobilized present in the reach BEFORE transfer
Qbi_dep{1} = zeros([size(Network.II),length(psi_3S)]);
Qbi_tr{2} = Qbi_tr{1};
Qbi_dep{2} = Qbi_dep{1};

Q_out{1} = zeros( length(Network.II), length(psi_3S));

% Routing scheme  

for t = 2:timescale-1
    
    %variables initialization
    Qbi_tr{2} = zeros( size(Qbi_tr{1}) );
    Qbi_dep{2} = zeros( size(Qbi_dep{1}) );
    Qbi_mob{1} = zeros( size(Qbi_mob{1}) );
    Q_out{t} = zeros( size(Q_out{1}) );
        
    %calculate new water dept for all reaches
    %Manning
    h = (Q(t,:).*n_Man./(Wac.*sqrt( Slope ))).^(3/5);
    v = 1./n_Man.*h.^(2/3).*sqrt( Slope ); 
  
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
        tr_cap = Engelund_Hansen_tr_cap(Fi_r_reach' , Slope(n) , Wac(n), v(n) , h(n) ,psi_3S)' .* 24.*60.*60;

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

    [ v_sed ] = velocity_EH( Slope , Wac , v , h , minvel, phi ,psi_3S) ;

    for n = NH

        V_mob = squeeze(Qbi_mob{1}(:,n,:));

        if sum(V_mob,'all') >0

            [reach_dest, setout] = sed_transfer_fast( n , v_sed.*(60*60*24) , Lngt, Network );
            
            % i sum the volumes transported in reach n with all the other
            % volumes mobilized by all the other reaches in time t
            
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

end