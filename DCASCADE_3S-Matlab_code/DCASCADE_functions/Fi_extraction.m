function [ Fi_r ] = Fi_extraction(  D50  , s ,psi )
%FI_EXTRACTION returns the sediment distribution for the considered sediment
%classes based on the D50 and the inverse value of spread s
%
%The function uses the Rosin sediment distribution (Shih and Komar, 1990)
%
%INPUT : 
% D50   = vector containing the D50 value for all reaches of the network
% psi   = vector containing the mean grain size of each sediment class in phi, it must be coherent with the sediment distribution in Fi_Sn (same number of sediment class)  
% s     = inverse measure of the spread of the GSD curve
%
%OUTPUT :
% Fi_Sn = Frequency of sediment for each sediment class, defined for each source node. 

%% GSD calculation
% convert D50 in mm
D50_finer = D50 * 1000;

%sediment classes diameter (mm)
dmi = 2.^(-psi); 

% find k parameter
k = D50_finer./((-log(1-50/100))^(1/s));

% find Fi_r
F = 1 - exp(-(flip(dmi)./k).^s);
F(: , size(F,2)  ) = 1;

Fi_r  = flip([F(:,1),diff(F,1,2)],2);


end