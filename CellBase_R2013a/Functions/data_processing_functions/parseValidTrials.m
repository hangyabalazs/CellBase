function valid_trials2 = parseValidTrials(VE,event,valid_trials)
%PARSEVALIDTRIALS   Deal with different forms of 'valid_trials'.
%   PARSEVALIDTRIALS is a helper function for input argument handling.
%   VALID_TRIALS = PARSEVALIDTRIALS(VE,EVENT,VALID_TRIALS) checks for NaNs
%   in the field corresponding to the EVENT in the events structure (VE,
%   stimulus or trial events). It handles VALID_TRIALS input in both
%   logical and index set format.
%
%   See also FILTERTRIALS.

% Remove NaNs from 'valid_trials'
if isequal(valid_trials,'all')
    valid_trials2 = ~isnan(VE.(event));
else
    if islogical(valid_trials)
        valid_trials2 = valid_trials & ~isnan(VE.(event));
    else
        valid_trials2 = intersect(valid_trials,find(~isnan(VE.(event))));
    end
end