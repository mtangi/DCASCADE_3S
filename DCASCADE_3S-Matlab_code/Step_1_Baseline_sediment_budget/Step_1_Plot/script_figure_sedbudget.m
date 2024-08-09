%% figure sed yield sensitivity analysis

load('DCASCADE_3S-Matlab_code/Step_1_Baseline_sediment_budget/Step_1_Result/output_baseline_sed_budget.mat')

%% plot figure
Fsize = 20;

sed_yield = 15:5:60;
sed_D50 = [1 0.75 0.5 0.25 0.1 0.075 0.05 0.025 ];
sed_partID = [1 2];

colrange = flip(jet(length(sed_D50)));

txt = cellstr(num2str([1:length(data_output_allcomb)]')); %txt = n days of flushing

accepted_range_y = [17 25]; %accepted range of sed.yield according to field data

figure
%plot acceptable range       
rectangle('Position',[accepted_range_y(1) 0  diff(accepted_range_y) 100],'FaceColor',[0.8 0.98 0.98] ,'EdgeColor',[0.7 0.88 0.88],'LineStyle' ,'--','LineWidth',3);
hold on 

j=1;

% fake point for legend

h = area([-1 -1],[0.1 0.1],'LineStyle','none');
h(1).FaceColor = [0.8 0.98 0.98];

plot(-1,-1,'Color','none');

sz = max(sed_yield)/5;

plot(-1,-1,'Color','none'); %empty spot in legend

for i=1:2:length(sed_yield)
    sz = sed_yield(i)/5;
    plot(-1,-1,'ro','MarkerSize',sz ,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    hold on
end

% plot figure
for i=1:length(data_output_allcomb(:,2 ))
    
    if ~isempty(data_output_allcomb{i,2}) && data_output_allcomb{i, 2}(3)==1 && data_output_allcomb{i, 2}(2)>0.01
        
        data(i,1) = mean(data_output_allcomb{i, 1}{2, 2}(2:end)); %mean yield
          
        data(i,2) = mean(data_output_allcomb{i, 1}{3, 2}(2:end)); %mean D50
        
        Fi_mean = mean(data_output_allcomb{i, 1}{4, 2}(:,2:end),2); % sand percentage
        data(i,2) = sum(Fi_mean(1:3))*100;
        
        if sum(Fi_mean(1:3)) > 0.1 && data(i,1) >= accepted_range_y(1) && data(i,1) <= accepted_range_y(2)
            solution_ok(j) = i;
            j=j+1;
        end
    
        sz = data_output_allcomb{i, 2}(1)*3; %size for sed. yield
        
        col = colrange(data_output_allcomb{i, 2}(2)==sed_D50,:); % color for D50in
        
        if data_output_allcomb{i, 2}(3)==2 %marker for sed. delivery method
            mrk = 'd'; else mrk = 'o';
        end
        
            scatter(data(i,1),data(i,2),sz,mrk,'filled','MarkerFaceColor',col,'MarkerEdgeColor','k')
        
        hold on
        
    end
      
end

txtlg = [ 'Acceptable range' ;' ' ;'Input sed. yield [Mt/y]'; cellstr(num2str(sed_yield(1:2:end)'))]; %txt = n days of flushing

legend(txtlg)

cmap = colormap(flip(jet(length(sed_D50))) ) ; %Create Colormap
cbh = colorbar ; %

caxis( [ min(sed_D50) max(sed_D50)] )

cbh.Ticks = linspace( min(sed_D50), max(sed_D50) , length(sed_D50)+1 )+0.05 ; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(sed_D50) ;

xlim([0 35])
ylim([0 100])
ylabel(cbh, 'Input D50 [mm]') ;

clear colortype c  cb ic nameclb varplot objplot order varid

grid on
grid minor
ylabel('d(3S)[mm]')
ylabel(' Sand and gravel content in the sediment delivered [%]')

xlabel('Sediment delivery to the Mekong [Mt/yr]')

set(gca,'FontSize',Fsize)
set(gcf,'color','w');

%% plot figure with pie chart for example scenarios
Fsize = 20;

examplepieID = [6 26 46 66 76];
    
sed_yield = 15:5:60;
sed_D50 = [1 0.75 0.5 0.25 0.1 0.075 0.05 0.025 ];
sed_partID = [1 2];

colrange = flip(jet(length(sed_D50)));

txt = cellstr(num2str([1:length(data_output_allcomb)]')); %txt = n days of flushing

accepted_range_y = [17 25]; %accepted range of sed.yield according to field data

figure
sb = subplot(4,5,[1:10]);
pos = get( sb, 'Position' );
pos(2) = 0.58;
set( sb, 'Position', pos ) ;

%plot acceptable range       
rectangle('Position',[accepted_range_y(1) 0  diff(accepted_range_y) 100],'FaceColor',[0.8 0.98 0.98] ,'EdgeColor',[0.7 0.88 0.88],'LineStyle' ,'--','LineWidth',3);
hold on 

j=1;

% fake point for legend

h = area([-1 -1],[0.1 0.1],'LineStyle','none');
h(1).FaceColor = [0.8 0.98 0.98];

plot(-1,-1,'Color','none');

sz = max(sed_yield)/5;

plot(-1,-1,'Color','none'); %empty spot in legend

for i=1:2:length(sed_yield)
    sz = sed_yield(i)/5;
    plot(-1,-1,'ro','MarkerSize',sz ,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    hold on
end

% plot figure
for i=1:length(data_output_allcomb(:,2 ))
    
    if ~isempty(data_output_allcomb{i,2}) && data_output_allcomb{i, 2}(3)==1 && data_output_allcomb{i, 2}(2)>0.01
        
        data(i,1) = mean(data_output_allcomb{i, 1}{2, 2}(2:end)); %mena yield
          
        data(i,2) = mean(data_output_allcomb{i, 1}{3, 2}(2:end)); %mean D50
        
        Fi_mean = mean(data_output_allcomb{i, 1}{4, 2}(:,2:end),2); % sand percentage
        data(i,2) = sum(Fi_mean(1:3))*100;
        
        if sum(Fi_mean(1:3)) > 0.1 && data(i,1) >= accepted_range_y(1) && data(i,1) <= accepted_range_y(2)
            solution_ok(j) = i;
            j=j+1;
        end
    
        sz = data_output_allcomb{i, 2}(1)*3; %size for sed. yield
        
        col = colrange(data_output_allcomb{i, 2}(2)==sed_D50,:); % color for D50in
        
        if data_output_allcomb{i, 2}(3)==2 %marker for sed. delivery method
            mrk = 'd'; else mrk = 'o';
        end
        
        scatter(data(i,1),data(i,2),sz,mrk,'filled','MarkerFaceColor',col,'MarkerEdgeColor','k')
        dx = 0.001; dy = 0.001; % displacement so the text does not overlay the data points
        %text( data(i,1)+dx, data(i,2)+dy, txt(i) );
        hold on

        hold on
        
    end
      
end

txtlg = [ 'Acceptable range' ;' ' ;'Input sed. yield [Mt/y]'; cellstr(num2str(sed_yield(1:2:end)'))]; %txt = n days of flushing

legend(txtlg)

cmap = colormap(flip(jet(length(sed_D50))) ) ; %Create Colormap
cbh = colorbar ; %

caxis( [ min(sed_D50) max(sed_D50)] )

cbh.Ticks = linspace( min(sed_D50), max(sed_D50) , length(sed_D50)+1 )+0.05 ; %Create 8 ticks from zero to 1
cbh.TickLabels = num2cell(sed_D50) ;

xlim([0 35])
ylim([0 100])
ylabel(cbh, 'Input D50 [mm]') ;

clear colortype c  cb ic nameclb varplot objplot order varid

grid on
grid minor
ylabel('d(3S)[mm]')
ylabel(' Sand and gravel content [%\Theta(3S)]')

xlabel('Sediment delivery to the Mekong (Mt/y)')

set(gca,'FontSize',Fsize)
set(gcf,'color','w');

% plot pie charts

Fi_in = cell2mat(cellfun(@(x)x{5,2},data_output_allcomb(examplepieID,1),'UniformOutput',0))';
Fi_out = cell2mat(cellfun(@(x)mean(x{4,2},2),data_output_allcomb(examplepieID,1),'UniformOutput',0)');

piefontsize = 15;

for i=1:size(Fi_in,2)
    s = subplot(4,5,i+10);
    
    p = pie(Fi_in(:,i), []);

    cmap = colormap(s, parula) ; %Create Colormap
    
    title(['D50_{in} = ' num2str( data_output_allcomb{examplepieID(i),2}(2) ) ' mm'],'FontSize',piefontsize+2)

    %title(['D50_{in} = ' num2str( sed_input_param_comb{examplepieID(i),2}(2) ) ' mm' ' - ' 'Q_{s-in} = ' num2str( sed_input_param_comb{examplepieID(i),2}(1) ) ' Mt/y' ],'FontSize',piefontsize+2)

    set(findobj(p,'type','text'),'FontSize',piefontsize)
end

for i=1:size(Fi_in,2)
    s=subplot(4,5,i+15);
    p = pie(Fi_out(:,i), []);

    cmap = colormap(s, parula) ; %Create Colormap

    set(findobj(p,'type','text'),'FontSize',piefontsize)
    %title({'';' '})
    
    if i==3
      legend({'Gravel', 'C. Sand ','F. Sand ' 'C. Silt ', 'F. Silt '},'Location','southoutside','Orientation','horizontal','FontSize',Fsize)      
    end
end

set(gcf,'color','w');
set(gca,'FontSize',15)

%clear data_D50 data_yield data_D50_in p Fi_in Fi_out

