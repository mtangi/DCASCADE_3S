function [st_m,m_y] = findFLmonth(variables)

% find the actual month when flushing is performed according to the
% startmonth and minflow parameters

%%
load('Q_data_3Sborg.mat')
load('DamDatabase_3S.mat')

st_m = zeros(size(variables(:,2)));

for s=1:size(variables,1)
    
    
    startmonth = variables(s,2);
    minflow = variables(s,3).* max(Q(:,[DamDatabase_active.reach_UP]),[],1) ;
    
    
    for i=1:length(Q)/365
                
        for d=1:length([DamDatabase_active.reach_UP])
            
            if ~isempty(find(Q(i+(365)*(i-1):365*i-1,[DamDatabase_active(d).reach_UP]) > minflow(d),1))
                m_y{s}(i,d) = dates_Q(2,find(Q(i+(365)*(i-1):365*i-1,[DamDatabase_active(d).reach_UP]) > minflow(d),1));
            else
                m_y{s}(i,d) = 1;
            end
        end
    end
    
    st_m(s) = max(round(mean(m_y{s},'all')),startmonth);
    
end


end

