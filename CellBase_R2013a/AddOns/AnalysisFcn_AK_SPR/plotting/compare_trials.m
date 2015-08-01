function [values, values_se, accu, num_trials] = compare_trials(EpochRate,TE,Analysis,COMP,valid_trials);
%
% COMPARE_TRIALS
%
%
% AK 07/1
 


%if valid_trials doesnt contain trial number then convert
if isbinary(valid_trials) 
     valid_trials = find(valid_trials);
end 

if iscell(COMP)
    NumConditions = length(COMP);
    for iC = 1:NumConditions
        trial_pos{iC} = intersect( COMP{iC}, valid_trials);
    end %i
else % matrix
   'problem';
end

accu = EpochRate*NaN;

if ~isempty(Analysis) 
    to_evaluate =  [strrep(Analysis,'trials','trial_pos{iC}') ';'];
    NumValues = 2;
else
    NumValues = 1;
end


values    = nan(NumConditions,NumValues);
values_se    = nan(NumConditions,NumValues);

for iC = 1:NumConditions
    
    NUM_TRIALS = length(trial_pos{iC});  % doesn't account for variable windows
    
    if ~isempty(trial_pos{iC})
        values(iC,1)     = nanmean(EpochRate(trial_pos{iC}));
        if NUM_TRIALS > 1
            values_se(iC,1)  = nanstd(EpochRate(trial_pos{iC}))./sqrt(NUM_TRIALS-1);
        else
            values_se(iC,1) = NaN;
        end
        num_trials(iC)   = NUM_TRIALS;
        if NumValues > 1
            junk = eval(to_evaluate);
            accu(trial_pos{iC}) = nanmean(junk);
            values(iC,2)        = nanmean(junk);
            values_se(iC,2)     = nanstd(junk)./sqrt(NUM_TRIALS-1);
        end
    end
end %iCOND