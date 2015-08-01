function [psth, spsth, spsth_se] = binraster2psth(binraster,dt,sigma,COMP,valid_trials)
%BINRASTER2PSTH    PSTH from binned spike raster.
%   [PSTH, SPSTH, SPSTH_SE] = BINRASTER2PSTH(BINRASTER,DT,SIGMA,COMP,VALID_TRIALS)
%   calculates PSTH, smoothed PSTH (SPSTH) and standard error of smoothed
%   PSTH (SPSTH_SE) from BINRASTER at time resolution DT. The output is
%   restricted to VALID_TRIALS (1 x NUM_TRIALS). COMP: M x NUM_TRIALS
%   matrix of integers (0...N_CONDITION), the matrix for partitioning
%   trials. SIGMA determines the smoothing kernel for the smoothed PSTH
%   (see SMOOTHED_PSTH).
%
%   See also SMOOTHED_PSTH, STIMES2BINRASTER and VIEWCELL2B.

%   Edit log: AK 06/1, AK 07/1, BH 6/23/11, 7/5/12

% If valid_trials doesn't contain trial number then convert
if isbinary(valid_trials) 
     valid_trials = find(valid_trials);
end

% Calculate positions
if iscell(COMP)
    nCOND = length(COMP);
    positions = cell(1,nCOND);
    for i = 1:nCOND
        positions{i} = intersect(COMP{i},valid_trials);
    end  % i
else   % matrix
    iCOND = 1;
    nCOND = size(COMP,1);
    for i = 1:nCOND
        CONDITIONS = unique(COMP(i,:));
        if length(CONDITIONS) > 2 && ~ismember(100,CONDITIONS)
            CONDITIONS = setdiff(CONDITIONS,0);   % conditions to compare, 0 doesn't count
        end
        for j = 1:length(CONDITIONS)
            positions{iCOND} = intersect(find(COMP(i,:)==CONDITIONS(j)),valid_trials); %#ok<AGROW>
            iCOND = iCOND + 1;
        end   % j
    end   % i
end    % if

% Preallocate output
psth = nan(nCOND,size(binraster,2));
psth_sd = psth;
spsth = psth;
spsth_se = psth;

% Calculate PSTH
for iCOND = 1:length(positions)
    NUM_TRIALS = length(positions{iCOND});  % doesn't account for variable windows
    if ~isempty(positions{iCOND})  
        psth(iCOND,:) = nanmean(binraster(positions{iCOND},:)) / dt; 
        psth_sd(iCOND,:) = nanstd(binraster(positions{iCOND},:)) / dt;
        spsth(iCOND,:) = smoothed_psth(psth(iCOND,:),dt,sigma);
        spsth_se(iCOND,:) = smoothed_psth(psth_sd(iCOND,:) ./ sqrt(NUM_TRIALS-1),dt,sigma);
    end
end   % iCOND