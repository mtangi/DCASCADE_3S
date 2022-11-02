%%
load('DCASCADE_3S-Matlab_code/Step_3_Reservoir_sediment_management/Part_1_flushing_design_optimization/Part_1_Result/output_Borg_result.mat')

objectives = output_Borg_result{1,2};
variables = output_Borg_result{1,3};

clear Jopt theta

%% round resulting parameters 

variables(:,1) = round(variables(:,1)); %duration
variables(:,2) = ceil(variables(:,2));  %startmonth

variables(:,4) = round(variables(:,4)); %annfreq
variables(variables(:,4)>3,4) = 4; %annfreq

variables(:,5) = round(variables(:,5)); %fllag

%% plot pareto front (fixed variables)

varidcl = 2;
varidsz = 4;
varidlv = 1;

% define id of variable to be shown as color 
nameclb = {'Average Flushing Month'};

%define type of market to be used and size
rangesz = [220 221];
mrkclass = {'d','s','<','o'};

% define id of variable to be shown as line width 
rangelv = [0.5 4.5];

%sort variable and objective according to varid

objplot = -objectives;
var_corrected = variables;
var_corrected(:,2) = findFLmonth(var_corrected); %change the starting month according to the average starting month for each dam for eac timestep

%define color for each portfolio

colortype = hsv(12);
colortype = abs(colortype - [1 1 1 ]); %negative of the color
c = colortype(var_corrected(:,varidcl),:);

%define size for each plot
sz = ((max(var_corrected(:,varidsz))+1 - var_corrected(:,varidsz))-min(var_corrected(:,varidsz))) / (max(var_corrected(:,varidsz))- min(var_corrected(:,varidsz))).*diff(rangesz) + min(rangesz);

mrktype = mrkclass(var_corrected(:,varidsz));

%define marker line wirdth for each plot
lv = (var_corrected(:,varidlv)-min(var_corrected(:,varidlv))) / (max(var_corrected(:,varidlv))- min(var_corrected(:,varidlv))).*diff(rangelv) + min(rangelv);

% define text to be plotted 
txt = cellstr(num2str([1:size(var_corrected,1)]')); %txt = solution ID
txt = cellstr(num2str(var_corrected(:,1))); %txt = n days of flushing


figure


% fake point for legend (frequency)
hold on
plot(-1,-1,'rd','MarkerSize',30 ,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])

plot(-1,-1,'rs','MarkerSize',30 ,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])

plot(-1,-1,'r<','MarkerSize',30 ,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])

%plot pareto front
for i=1:size(var_corrected,1)
    
    %scatter(objplot(i,1)/1E3,objplot(i,2),sz(i),c(i,:),'filled','LineWidth',lv(i),'Marker',mrktype{ i},'MarkerEdgeColor',[0.5 0.5 0.5]);
    scatter(objplot(i,1)/1E3,objplot(i,2),sz(i),c(i,:),'filled','LineWidth',0.5,'Marker',mrktype{ i},'MarkerEdgeColor',[0 0 0 ]);

    hold on
    %text( objplot(i,1)/1E3, objplot(i,2), txtmrk(i),'FontSize' ,15,'HorizontalAlignment','center' );
    
%     hold on
%     dx = 0.001; dy = 0.1; % displacement so the text does not overlay the data points
%     text( objplot(i,1)/1E3+dx, objplot(i,2)+dy, txt(i) );
 
end

% for i=unique(var_corrected(:,varidlv))'
% 
%     scatter(objplot(var_corrected(:,varidlv) == i,1)/1E3,objplot(var_corrected(:,varidlv) == i,2),sz(var_corrected(:,varidlv) == i),c(var_corrected(:,varidlv) == i,:),'filled','LineWidth',lv(i==unique(var_corrected(:,varidlv))'),'MarkerEdgeColor','k');
%     hold on
%     
%     dx = 0.001; dy = 0.1; % displacement so the text does not overlay the data points
%     text( objplot(var_corrected(:,varidlv) == i,1)/1E3+dx, objplot(var_corrected(:,varidlv) == i,2)+dy, txt(var_corrected(:,varidlv) == i) );
%     
% end

scatter( 1.39428137224075 , 13.4321757848378,1000,[0 0 0],'filled','Marker','hexagram','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth' , 0.5);

    
xlabel('Hydropower generation [10^3 GWh/yr]')
ylabel('Sediment load at network outlet [Mt/yr]')

%set colorbar

cmap = colormap(colortype ) ; %Create Colormap
cbh = colorbar ; %

caxis( [ 1 13 ] )

cbh.Ticks =1.5:13 ; %Create 12 ticks from 1 to 12
cbh.TickLabels = {'Gen','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',""} ;

ylabel(cbh, nameclb) ;
ylim([10,26]);
xlim([1.2 1.4])
set(gca,'FontSize',22)
set(gcf,'color','w');

lgd=legend({'1yr                      .','2yr', '3yr' },'Location','southwest');

lgd.FontSize = 20;
lgd.Title.String = 'Flushing frequency';

clear colortype c  cb ic nameclb varplot objplot order varid lv sz rangelv mrktype mrkclass 
clear rangesz varidcl varidlv varidsz varcolor cmap i cbh a dx dy txt txtrmrkclass txtmrk


