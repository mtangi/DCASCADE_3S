function [] = plot_dam_features(dam_output,dates_Q,varargin)
% PLOT_DAM_FEATURES shows the evolution through time of the reservoir
% features, including water level, water input, standard release and
% release through spillways.

%% read additional inputs 

p = inputParser;
addOptional(p,'plot2axis', 3 );
addOptional(p,'timeinterval',[]);
addOptional(p,'LineWidth',1.5 );
addOptional(p,'FontSize',12);
addOptional(p,'titles',[]);
addOptional(p,'showlegend', 'on');

parse(p,varargin{:})

LineWidth = p.Results.LineWidth;
TI = p.Results.timeinterval;
FontSize = p.Results.FontSize;
titles_subplot = p.Results.titles;
showlegend = p.Results.showlegend;

%% load dam features

figure

set(gcf,'color','w');
DamDatabase_active = dam_output{1,2};
ResVolume = dam_output{3,2};
release = dam_output{4,2};
sed_storage = dam_output{6,2};
Q =  dam_output{9,2};

n_dams = length(DamDatabase_active);

ts = size(ResVolume,1)-2;

%% set time intervale
if isempty(TI )
    TI = [2 ts];
end

%% plot dam features

dam_heigth = zeros(ts,size(release,2));

for d=1:n_dams
    sb = subplot(n_dams,1,d);
    pos = get(sb, 'position');

    %plot inflow 
    plot( datetime(dates_Q(:,1:ts)'),[Q(1:ts, DamDatabase_active(d).reach_UP)],'LineWidth',LineWidth)
    hold on

    %plot release (outflow)
    plot( datetime(dates_Q(:,1:ts)'),[release(1:ts,d)],'LineWidth',LineWidth)
    hold on

    %plot splillway release
    spillway = max(0,release(1:ts,d) - DamDatabase_active(d).design_discharge);
    plot( datetime(dates_Q(:,1:ts)'),spillway,'LineWidth',LineWidth)

    ylim([0, max([release;Q(1:ts, [DamDatabase_active.Node_ID])] ,[],'all')]);
    ylabel('Discharge [m^3/s]')

    yyaxis right
    
    % plot reservoir level
    dam_heigth(:,d) = Reservoir_LSWConversion_FL( [ResVolume(1:ts,d)] ,DamDatabase_active(d).WL_table, 3, 1);
 
    ylim([0, max(dam_heigth,[],'all')+0.5]);

    plot( datetime(dates_Q(:,1:ts)'), dam_heigth(:,d) ,'LineWidth',LineWidth)

    %add new title if specified
    if isempty(titles_subplot) 
        title([DamDatabase_active(d).Name ' - reach ' num2str(DamDatabase_active(d).Node_ID)])
    else      
        title( [titles_subplot{d} ] );
    end

    ylabel('Water level [m]')
    set(gca,'FontSize',FontSize)

    xlim( [ datetime(dates_Q(:,TI(1))')  datetime(dates_Q(:,TI(2))') ] )
    grid on

    set(gca, 'YGrid', 'off', 'XGrid', 'on')
end

if strcmp(showlegend,'on')
    legend({'Inflow [m^3/s]','Outflow [m^3/s]','Spillway [m^3/s]','Reservoir heigth [m]'},'Location','south','Orientation','horizontal','FontSize',FontSize)
end


set(gcf,'color','w');

end

