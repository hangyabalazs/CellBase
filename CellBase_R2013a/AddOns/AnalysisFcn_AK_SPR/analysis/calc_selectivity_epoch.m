function out = calc_selectivity_epoch(cellid,EpochName1,EpochName2,nboot,eventtype)
% 
% CALC_SELECTIVITY_EPOCH
%
% Calculates the area under the ROC measure for the differences in firing rate
% across two different epochs and the bootsrapped P-value.
% eventtype = 'stim' or 'behav'
% EpochName1 and EpochName2 should be used to decide eventtype in the
% future
%
% calc_selectivity_epoch(cellid,EpochName1,EpochName2,nboot)


% AK 10/05 SPR 2010-08-04

if nargin<4,
    eventtype='behav';
end
switch eventtype
    case 'stim'
        try
            ST = loadcb(cellid,'STIMSPIKES');
            TE = loadcb(cellid,'StimEvents');
        catch
            out=[NaN NaN NaN NaN NaN];
        end
    case 'behav'
        try
            ST = loadcb(cellid,'EVENTSPIKES');
            TE = loadcb(cellid,'TrialEvents');
        catch
            out=[NaN NaN NaN NaN NaN];
        end
end

if ~exist('out','var'),
    epoch_pos1 = findcellstr(ST.epochs(:,1),EpochName1);
    epoch_pos2 = findcellstr(ST.epochs(:,1),EpochName2);
    
    if (epoch_pos1 == 0) | (epoch_pos2 == 0)
        error('Epoch name not found');
    end
    
    RATE1 = ST.epoch_rates{epoch_pos1};
    RATE2 = ST.epoch_rates{epoch_pos2};
    
    try
        ValidTrials = selecttrial(TE,'OdorPokeValid & WaterPokeValid ');
    catch
        ValidTrials=1:length(TE.TrialStart);
    end
    
    [D,  P]    =  rocarea(RATE1(ValidTrials),RATE2(ValidTrials),'boot',nboot,'scale');
    
    out = [D P length(ValidTrials) nanmean(RATE1) nanmean(RATE2)];
else
end
 