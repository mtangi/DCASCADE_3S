function [output] = Reservoir_LSWConversion(input ,input_table, id_input, id_output)
%Reservoir_LSWConversion interpolates the current value of water level(L),
%surface(S) or volume (V) of a reservoir (according to id_output), given some
%standard values of these parameters contained in table and an input value
%of a different parameter (defined by id_output).

% id_input and output
% 1: Supply level
% 2: Reservoir area
% 3: Reservoir volume

%% extract table values 
output = zeros(size(input));

for i=1:length(input)
column_input = input_table(:,id_input);
column_output = input_table(:,id_output);

pos_high = find(input(i)<column_input,1);
pos_low = pos_high - 1;

if sum(input(i)>column_input) >= length(column_input)
    warning('the input value is higher then the tabulated values');
end

%extract the output value by performing a linear interpolation between
%the lower and upper value of the input value, given by the table.

output(i) = column_output(pos_low)+( (input(i) - column_input(pos_low)) / (column_input(pos_high) - column_input(pos_low)) * (column_output(pos_high) - column_output(pos_low)) );

end

end

