function [Qbi_input] = initialize_sed_input(ReachData, Fi_yield, Q, dates_Q , psi, roundpar )
%%initialize_sed_input calculate the daily sed input in each reach, given
%%the annual sed.input.
% Annual sed.input is partitioned proportionally to the daily discharge.

%% variable initialization

clear Qbi_input; Qbi_input  = cell(size(Q,1),1); [Qbi_input{:}] = deal(zeros(length(ReachData), length(psi)));

dates_Q = dates_Q(:,1:min(length(Q),length(dates_Q))); %check if dates_Q and Q are of the same length
Q = Q(1:min(length(Q),length(dates_Q)),:);

years_ts = unique(dates_Q(1,:)); %years in the simulation

%sed yield in t km-2 y-1
rho_s = 2.650; %sediment density (t/m3)
sed_yield_year = [ReachData.sed_yield]' ./ rho_s .*[ReachData.directAd]'; %sediment yield (m3/*yr)

if isempty(roundpar)
    roundpar = rp;
end

if isempty(psi)
    roundpar = PSI;
end

%% calculate daily sed. yield from catchment in Qbi_input

sed_yield_day = zeros(length(ReachData),length(Q),length(psi));
j=1;

%for each year...
for i=1:length(years_ts)

    length_year = sum(dates_Q(1,:) == years_ts(i));

    %for each reach...
    for n=1:length(ReachData)
        daily_Q_perc_n = Q(dates_Q(1,:) == years_ts(i),n) ./ sum(Q(dates_Q(1,:) == years_ts(i),n)); %discharge for each day of the year [Kg/m3]
        sed_yield_year_n_class = sed_yield_year(n).* Fi_yield(n,:);   %annual sediment yield for reach n for each sed.class (m3/yr)

        sed_yield_day_tot_n = round(daily_Q_perc_n.*sed_yield_year(n).* Fi_yield(n,:),roundpar); %daily sediment yield for reach n for each sed.class (m3/day)

        %If i just round the sediment input, I risk losing some sediment volumes,
        %espacially for less represented sediment classes.
        %Thus, i need to add the sediment volume lost in the round
        %operation back to the inputs. I do so by adding the difference
        %(diff) between the rounded and non-rounded yield back to the
        %rounded yield matrix. 
        diff = round(sed_yield_year_n_class - sum(sed_yield_day_tot_n));

        %for each sed.class...
        for d=1:length(diff)
            [~,index] = sortrows(sed_yield_day_tot_n,d,'desc'); % i select the days in the year with the highest yield
            sed_yield_day_tot_n(index(1:abs(diff(d))),d) = sed_yield_day_tot_n(index(1:abs(diff(d))),d)+diff(d)/abs(diff(d)); % I add or subtract the difference in these days, spread so that the added or sub. value is 1m3 per day
        end

        if any(sed_yield_day_tot_n<0,'all')
            warning(['Negative value in reach ' num2str(n) ': ' num2str(min(sed_yield_day_tot_n,[],'all')) ]) 
            sed_yield_day_tot_n(sed_yield_day_tot_n<0) = 0;
        end

        sed_yield_day(n,j:j+length_year-1,:) = sed_yield_day_tot_n;

    end

    j=j+length_year;
end

Qbi_input = cell(length(Q),1);

for i=1:length(Qbi_input)
    Qbi_input{i} = squeeze(sed_yield_day(:,i,:));   
end
