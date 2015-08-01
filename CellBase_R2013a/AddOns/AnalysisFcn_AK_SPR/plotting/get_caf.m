function [caf, caf_se, xbins, nt] = get_caf(RATE,TE,Partitions,valid_trials,NumBins,MinPoints,Method);
%
% GET_CAF
%

if nargin < 6
    MinPoints = 5;
    Method = 'fixed_bins';
end

%Partitions = {'Correct','Error'};
CT = partition_trials(TE,Partitions);
posCORR = intersect(CT{1},valid_trials);
posERR  = intersect(CT{2},valid_trials);

[caf, xbins, nt] = histratio(RATE(posCORR),RATE(posERR),NumBins,Method,'minpoints',MinPoints);
caf_se = sqrt(caf-caf.^2) ./ sqrt(nt-1);  %binomiaml errorbar
