function values = get_tuning(TrialData,Analysis,COMP,valid_trials);


% --> exectute expression on subset of trials
%   length(unique(TE.OdorRatio))
%   mean(TE.Correct(valid_trials))
%   min(TE.WaterWaitDur(valid_trials))
%   
% --> execture expression on epoch_rates
%     mean(epoch_rates(my_epoch))
%     
% --> execute expression on spikes
%     TriggerEvent
%     window
%     expression
%      [values, values_se, accu] = compare_trials(EpochRate,TE,Analysis,COMP,valid_trials);
% %

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
   error('get_tuning: problem');
end


NumAnalysis = length(Analysis);
values = nan(NumConditions,NumAnalysis);

for iC = 1:NumConditions
    
    NumTrials = length(trial_pos{iC});  % doesn't account for variable windows
    
    data = TrialData(trial_pos{iC});
    
    if ~isempty(data)
        for iA = 1:NumAnalysis
            values(iC,iA)  = eval(Analysis{iA});
        end
    end
end %iCOND