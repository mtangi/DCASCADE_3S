load('DCASCADE_3S-Matlab_code/Step_2_Reservoir_siting/Step_2_Result/output_damimpact.mat')

sed_storage_loss=cell(1,3);

for i=1:length(output_damimpact)
    
    scenario = output_damimpact{i,3}(3);
    
    if scenario ~=0
        
        if isempty(sed_storage_loss{scenario})
            sed_storage_loss{scenario} =  1- output_damimpact{i, 2}{5, 2}(end-3,:) ./ output_damimpact{i, 2}{5, 2}(1,:);
        else
            sed_storage_loss{scenario} =  [sed_storage_loss{scenario} ; 1- output_damimpact{i, 2}{5, 2}(end-3,:) ./ output_damimpact{i, 2}{5, 2}(1,:)];
        end
    end
    
end

sed_storage_loss_full = cell2mat(sed_storage_loss).*100;

%% boxplot

figure

boxplot(sed_storage_loss_full,'Labels',{'LSS2 (LSS2Only)','LSS2 (FullDam)', 'LSP3' , 'LSS3' , 'US1','LSS2-II','LSP2'})

ylabel('Reservoir storage loss after 10 yrs [%] ')
xlabel('Dam')

color_scenarios = parula(3);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    switch j
        case 1
            colors = color_scenarios(1,:);
        case 4
            colors = color_scenarios(2,:);
        case 7
            colors = color_scenarios(3,:);
    end
    patch(get(h(j),'XData'),get(h(j),'YData'),colors,'FaceAlpha',.5);
end

h = get(gca,'Children');
set(gca,'Children',[h(8) h(5) h(7) h(3) h(4) h(6) h(2) h(1)])
legend({'LSS2Only','FullDam','FullDam\_ alt'},'Location','northwest')

ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

set(gcf,'color','w');
set(gca,'FontSize',15)

