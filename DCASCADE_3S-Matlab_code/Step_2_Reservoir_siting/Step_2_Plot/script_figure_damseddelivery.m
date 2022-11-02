%% prepare data for bar plot
load('DCASCADE_3S-Matlab_code/Step_2_Reservoir_siting/Step_2_Result/output_damimpact.mat')
 
load('network_data_3Sborg.mat')

river_reach{1,1} = 'Se Kong';
river_reach{1,2} = [find(Network.Upstream.Distance{74} ~= Inf),75];
river_reach{2,1} = 'Se San';
river_reach{2,2} = [find(Network.Upstream.Distance{162} ~= Inf)];
river_reach{3,1} = 'Sre Pok';
river_reach{3,2} = [find(Network.Upstream.Distance{243} ~= Inf)];
river_reach{4,1} = 'Se San - Sre Pok confluence';
river_reach{4,2} = [163:167];

%id of the dam scenario 
dam_scenario_ID = cell2mat(cellfun(@(x)x(3),output_damimpact(:,3),'UniformOutput',0));

% find data with no dams

for i=0:length(unique(dam_scenario_ID))-1

    data_outcum_full = permute(reshape(cell2mat(cellfun(@(x)x{1, 2},output_damimpact(dam_scenario_ID==i,1),'UniformOutput',0))' , 5, 5, sum(dam_scenario_ID==i) ),[2 1 3]) ;
    data_outcum = data_outcum_full;
    data_outcum(:,2,:) = sum(data_outcum_full(:,1:3,:),2);
    data_outcum(:,3,:) = sum(data_outcum_full(:,4:5,:),2);
    data_outcum(:,1,:) = sum(data_outcum_full,2);

    data_outcum(:,4:5,:) = [];
    
    outcum_provenance{i+1,1} = round(mean(data_outcum,3))*2.6./1000000;
    outcum_provenance{i+1,2} = data_outcum*2.6./1000000;
end

%% bar plot of sediment provenance with difference in sed.size (with uncertanty bar) as waterfall

allfontsize = 22;
piefontsize = 13;

data = outcum_provenance(:,1);

figure
set(gcf,'color','w');

titles = {'NoDams', 'LSS2Only','FullDam','FullDam\_alt'};

colorscale = parula(3);
% colorscale(1,:) = [0.2021    0.4788    0.9911];
% colorscale(2,:) = [0.1540    0.5902    0.9218];
X = categorical(river_reach(1:3,1)');
X = reordercats(X,river_reach(1:3,1)');
    

for i=1:length(data)

    a= i-1;
    coord = [a*3+1+a:a*3+3+a a*3+16+a:a*3+18+a a*3+31+a:a*3+33+a];
    subplot(3,15,coord)

    model_series = data{i}(1:3,:); 
    model_series(:,3) = model_series(:,1); 
    clear model_error
    for c=1:size(data{i},2)
        %err(:,c) = prctile(squeeze( sum(outcum_provenance{i,2}(1:3,c,:),2) )'*2.6./1000000,95) - prctile( squeeze( sum(outcum_provenance{i,2}(1:3,c,:),2) )' *2.6./1000000,5);
        model_error(:,c) = (prctile(squeeze( sum(outcum_provenance{i,2}(1:3,c,:),2) ),95,2)  - prctile(squeeze( sum(outcum_provenance{i,2}(1:3,c,:),2) ),5,2) );
        
        if c==size(data{i},2); a=1; else; a=c; end
        error_middle(:,c) = (prctile(squeeze( sum(outcum_provenance{i,2}(1:3,c,:),2) ),95,2)  - prctile(squeeze( sum(outcum_provenance{i,2}(1:3,c,:),2) ),5,2)  )/2 + prctile(squeeze( sum(outcum_provenance{i,2}(1:3,a,:),2) ),5,2);
    end

    b = bar(model_series, 'grouped','EdgeColor','k');

        title(titles(i))

    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(model_series);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for j = 1:nbars
        x(j,:) = b(j).XEndPoints;
    end
    x=x';
    % Plot the errorbars
    errorbar(x,error_middle,model_error/2,'k','linestyle','none');

    %plot bar to cancel to make bar "float"
    zero_bar = zeros(size(model_series));
    zero_bar(:,3) = model_series(:,2);
    
    b = bar(zero_bar, 'grouped','EdgeColor','w','FaceColor','w');

    set(gca,'xticklabel',X)

    h=gca; h.XAxis.TickLength = [0 0];
    
    ylim([0 10.5])
    %xlabel('River')
    set(gca,'FontSize',allfontsize)

    if i==1
        
        ylabel('Sediment load at network outlet [Mt/yr]')
    end
    
    %plot line to connect bats
    plotlinex = x(:,2:3)';
    plotlinex(1,:) = plotlinex(1,:)-0.09;
    plotlinex(2,:) = plotlinex(2,:)+0.09;
    plot(plotlinex, [model_series(:,2),model_series(:,2)]','k' ,'LineWidth',1.2)
    

   hold off
     
end

legend({'Total','Gravel+Sand','Silt'},'Location','southoutside','Orientation','horizontal','FontSize',allfontsize)

clear X data outcum_provenance_river_old allfontsize piefontsize c h b a x

