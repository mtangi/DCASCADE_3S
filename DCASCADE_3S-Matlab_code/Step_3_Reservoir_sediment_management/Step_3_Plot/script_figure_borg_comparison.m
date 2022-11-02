%% load data on average scenario
load('DCASCADE_3S-Matlab_code/Step_3_Reservoir_sediment_management/Part_2_sensitivity_analysis/Part_2_Result/output_sensitivity.mat')

%% plot pareto front (with uncertanty bars)

sz = 150;
mrkclass = {'d','s','<','hexagram'};
varidsz = 4;
colortype = hsv(3);
mec = [0.5 0.5 0.5]; %'MarkerEdgeColor',
lwd = 0.5; %linewidth 

colortype_err = (1-hsv(3))*0.6 + hsv(3);
colortype_err = [1,0.600,0.600;0.500,0.9,0.50;0.600,0.600,1];

figure

% fake point for legend (equations)

% scatter(0,0,sz,colortype(1,:),'filled','Marker',mrkclass{ 1},'MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on
% scatter(0,0,sz,colortype(2,:),'filled','Marker',mrkclass{ 1},'MarkerEdgeColor',[0.5 0.5 0.5]);
% scatter(0,0,sz,colortype(3,:),'filled','Marker',mrkclass{ 1},'MarkerEdgeColor',[0.5 0.5 0.5]);
% scatter(0,0,sz,'k','filled','Marker',mrkclass{ 4},'MarkerEdgeColor',[0.5 0.5 0.5]);
 
% fake point for legend (frequency)
hold on
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor',colortype(1,:),'MarkerEdgeColor',[1 1 1 ])
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor',colortype(2,:),'MarkerEdgeColor',[1 1 1])
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor',colortype(3,:),'MarkerEdgeColor',[1 1 1])
plot(-1,-1,'hexagram','MarkerSize',30 ,'MarkerFaceColor','k','MarkerEdgeColor',[0.5 0.5 0.5])

%plot errors
for s=1:size(output_sensitivity,1)
    
    objplot = -output_sensitivity{s,2};
    var_corrected = output_sensitivity{s,3};
    
    %plot pareto front
    for i=1:size(output_sensitivity{s,4},1)-1
             
        findValues = cellfun(@(x)x(i+1,:),output_sensitivity(s,4),'UniformOutput',0);
        findValues = -cell2mat(findValues{1}');
        
        findValues(:,1) =  findValues(:,1)/1E3;
                 
        %plot error bar for generation
        maxvalues = prctile(findValues,90,1);
        minvalues = prctile(findValues,10,1);
                
        yneg = max( objplot(i,2) - minvalues(2) ,0);
        ypos = max( maxvalues(2) - objplot(i,2) ,0);
        xneg = max( objplot(i,1)/1E3 - minvalues(1) ,0);
        xpos = max( maxvalues(1) - objplot(i,1)/1E3 ,0);
        
        %plot error bar for generation
        
        if round(var_corrected(i,4)) >= 4
            color_err = [0 0 0];
        else
            color_err = colortype_err(s,:);
        end
                    
        errorbar(objplot(i,1)/1E3 , objplot(i,2),yneg,ypos,xneg,xpos,'.','CapSize' ,2,'LineWidth' ,0.5,'Color',color_err)

        hold on
         
    end


end


%plot central points
for s=1:size(output_sensitivity,1)
    
    objplot = -output_sensitivity{s,2};
    var_corrected = output_sensitivity{s,3};
    
    var_corrected(:,1) = round(var_corrected(:,1)); %duration
    var_corrected(:,2) = ceil(var_corrected(:,2));  %startmonth

    var_corrected(:,4) = round(var_corrected(:,4)); %annfreq
    var_corrected(var_corrected(:,4)>3,4) = 4; %annfreq

    var_corrected(:,5) = round(var_corrected(:,5)); %fllag

    var_corrected(:,2) = findFLmonth(var_corrected); %change the starting month according to the average starting month for each dam for eac timestep

    var_corrected(:,4) = round(var_corrected(:,4));
    
    mrktype = mrkclass(var_corrected(:,varidsz));

    %plot pareto front
    for i=1:size(output_sensitivity{s,4},1)-1
               
        if round(var_corrected(i,4)) >= 4
            color_scatter = [0 0 0];
        else
            color_scatter = colortype(s,:);
        end
        
        hold on
        
        if s==4
            sz=1000;
            mec = [ 1 1 1];
            lwd = 1.5;
        end
        scatter(objplot(i,1)/1E3,objplot(i,2),sz,color_scatter,'filled','Marker',mrktype{ i},'MarkerEdgeColor',mec,'LineWidth' , lwd);
        
%         text( objplot(i,1)/1E3, objplot(i,2), txtmrk(i),'FontSize' ,15,'HorizontalAlignment','center' );
% 
%         hold on
%         dx = 0.001; dy = 0.1; % displacement so the text does not overlay the data points
%         text( objplot(i,1)/1E3+dx, objplot(i,2)+dy, txt(i) );

    end

    hold off
end

xlabel('Hydropower generation [10^3 GWh/yr]')
ylabel('Sediment load at network outlet [Mt/yr]')

set(gca,'FontSize',22)

lgd=legend({'ATK','E&H', 'E&H - Width correction' , 'No flushing' },'Location','southwest');
lgd=legend({'Atkinson','Engelund', 'Engelund - Width correction' , 'No flushing' },'Location','southwest');


lgd.FontSize = 20;
lgd.Title.String = 'Flushing tr.cap. equations';

set(gcf,'color','w');

ylim([10 26])
xlim([1.2 1.4])

%% plot pareto front (with uncertanty bars and LSS2)

sz = [50 50 50 1000 1200];
lwd = [0.5 0.5 0.5 2 2];
mrkclass = {'d','s','<','hexagram','hexagram'};
varidsz = 4;
colortype = hsv(3);

colortype_err = (1-hsv(3))*0.6 + hsv(3);
colortype_err = [1,0.600,0.600;0.500,0.9,0.50;0.600,0.600,1];

figure
 
% fake point for legend (frequency)
hold on
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor',colortype(1,:),'MarkerEdgeColor',[0.5 0.5 0.5])
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor',colortype(2,:),'MarkerEdgeColor',[0.5 0.5 0.5])
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor',colortype(3,:),'MarkerEdgeColor',[0.5 0.5 0.5])
plot(-1,-1,'ro','MarkerSize',30 ,'MarkerFaceColor','k','MarkerEdgeColor',[0.5 0.5 0.5])
plot(-1,-1,'hexagram','MarkerSize',30 ,'MarkerFaceColor','k','MarkerEdgeColor',[0.5 0.5 0.5])

%plot errors
for s=1:3%size(Borg_results_sensitivity,1)
    
    objplot = -output_sensitivity{s,2};
    var_corrected = output_sensitivity{s,3};
    
    %plot pareto front
    for i=1:size(output_sensitivity{s,4},1)-1
             
        findValues = cellfun(@(x)x(i+1,:),output_sensitivity(s,4),'UniformOutput',0);
        findValues = -cell2mat(findValues{1}');
        
        findValues(:,1) =  findValues(:,1)/1E3;
                 
        %plot error bar for generation
        maxvalues = prctile(findValues,90,1);
        minvalues = prctile(findValues,10,1);
                
        yneg = max( objplot(i,2) - minvalues(2) ,0);
        ypos = max( maxvalues(2) - objplot(i,2) ,0);
        xneg = max( objplot(i,1)/1E3 - minvalues(1) ,0);
        xpos = max( maxvalues(1) - objplot(i,1)/1E3 ,0);
        
        %plot error bar for generation
        
        if round(var_corrected(i,4)) >= 4
            color_err = [0 0 0];
        else
            color_err = colortype_err(s,:);
        end
                    
        errorbar(objplot(i,1)/1E3 , objplot(i,2),yneg,ypos,xneg,xpos,'.','CapSize' ,2,'LineWidth' ,lwd(s),'Color',color_err)

        hold on
         
    end

        objplot = -output_sensitivity{s,2};
    var_corrected = output_sensitivity{s,3};
    
    var_corrected(:,1) = round(var_corrected(:,1)); %duration
    var_corrected(:,2) = ceil(var_corrected(:,2));  %startmonth

    var_corrected(:,4) = round(var_corrected(:,4)); %annfreq
    var_corrected(var_corrected(:,4)>3,4) = 4; %annfreq

    var_corrected(:,5) = round(var_corrected(:,5)); %fllag

    var_corrected(:,2) = findFLmonth(var_corrected); %change the starting month according to the average starting month for each dam for eac timestep

    var_corrected(:,4) = round(var_corrected(:,4));
    
    if s==5
         var_corrected(:,4) = 5;
    end
        
    mrktype = mrkclass(var_corrected(:,varidsz));

    %plot pareto front
    for i=1:size(output_sensitivity{s,4},1)-1
               
        if round(var_corrected(i,4)) >= 4
            color_scatter = [0 0 0];
        else
            color_scatter = colortype(s,:);
        end
        
        hold on
        scatter(objplot(i,1)/1E3,objplot(i,2),sz(s),color_scatter,'filled','Marker',mrktype{ i},'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',0.05);
        
%         text( objplot(i,1)/1E3, objplot(i,2), txtmrk(i),'FontSize' ,15,'HorizontalAlignment','center' );
% 
%         hold on
%         dx = 0.001; dy = 0.1; % displacement so the text does not overlay the data points
%         text( objplot(i,1)/1E3+dx, objplot(i,2)+dy, txt(i) );

    end


end

hold on

for s=3:size(output_sensitivity,1)
    
    objplot = -output_sensitivity{s,2};
    var_corrected = output_sensitivity{s,3};
    
    %plot pareto front
    for i=1:size(output_sensitivity{s,4},1)-1
             
        findValues = cellfun(@(x)x(i+1,:),output_sensitivity(s,4),'UniformOutput',0);
        findValues = -cell2mat(findValues{1}');
        
        findValues(:,1) =  findValues(:,1)/1E3;
                 
        %plot error bar for generation
        maxvalues = prctile(findValues,90,1);
        minvalues = prctile(findValues,10,1);
                
        yneg = max( objplot(i,2) - minvalues(2) ,0);
        ypos = max( maxvalues(2) - objplot(i,2) ,0);
        xneg = max( objplot(i,1)/1E3 - minvalues(1) ,0);
        xpos = max( maxvalues(1) - objplot(i,1)/1E3 ,0);
        
        %plot error bar for generation
        
        if round(var_corrected(i,4)) >= 4
            color_err = [0 0 0];
        else
            color_err = colortype_err(s,:);
        end
                    
        errorbar(objplot(i,1)/1E3 , objplot(i,2),yneg,ypos,xneg,xpos,'.','CapSize' ,2,'LineWidth' ,lwd(s),'Color',color_err)

        hold on
         
    end

        objplot = -output_sensitivity{s,2};
    var_corrected = output_sensitivity{s,3};
    
    var_corrected(:,1) = round(var_corrected(:,1)); %duration
    var_corrected(:,2) = ceil(var_corrected(:,2));  %startmonth

    var_corrected(:,4) = round(var_corrected(:,4)); %annfreq
    var_corrected(var_corrected(:,4)>3,4) = 4; %annfreq

    var_corrected(:,5) = round(var_corrected(:,5)); %fllag

    var_corrected(:,2) = findFLmonth(var_corrected); %change the starting month according to the average starting month for each dam for eac timestep

    var_corrected(:,4) = round(var_corrected(:,4));
    
    if s==5
         var_corrected(:,4) = 5;
    end
        
    mrktype = mrkclass(var_corrected(:,varidsz));

    %plot pareto front
    for i=1:size(output_sensitivity{s,4},1)-1
               
        if round(var_corrected(i,4)) >= 4
            color_scatter = [0 0 0];
        else
            color_scatter = colortype(s,:);
        end
        
        hold on
        scatter(objplot(i,1)/1E3,objplot(i,2),sz(s),color_scatter,'filled','Marker',mrktype{ i},'MarkerEdgeColor',[0.5 0.5 0.5]);
        
%         text( objplot(i,1)/1E3, objplot(i,2), txtmrk(i),'FontSize' ,15,'HorizontalAlignment','center' );
% 
%         hold on
%         dx = 0.001; dy = 0.1; % displacement so the text does not overlay the data points
%         text( objplot(i,1)/1E3+dx, objplot(i,2)+dy, txt(i) );

    end


end
hold off


xlabel('Hydropower generation [10^3 GWh/yr]')
ylabel('Sediment load at network outlet [Mt/yr]')

set(gca,'FontSize',28)

%legend({'ATK','E&H', 'E&H - Width correction' , 'FullDam - No flushing', 'LSS2 ' },'Location','southwest');

set(gcf,'color','w');

ylim([7 26])
xlim([1.2 2.1])


xticks([1.25:0.25: 2])
yticks([5:5:25])