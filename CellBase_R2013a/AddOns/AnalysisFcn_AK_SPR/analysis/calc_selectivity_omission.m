function out = calc_selectivity_omission(cellid,EpochName,nboot,eventtype,lookback);
% OUT = CALC_SELECTIVITY_OMISSION(CELLID,EPOCHNAME,NBOOT,EVENTTYPE,LOOKBACK)
%
% Calculates the area under the ROC measure for selectivity for reward omitted trials in 
% requested epoch and a bootsrapped P-value.
%
% lookback = 1; will sort the trials by PreviousRewardOmit (default)
% calc_selectivity_omission(cellid,EpochName,nboot)

%SPR 2010-12-31

if strcmpi(cellid,'default')
    out = [strcat({'Domi_','Pomi_'},{EpochName,EpochName})];
    return;
end

if nargin<4,
    eventtype='behav';
end

if nargin<5,
    lookback = 1;
end

switch eventtype
    case 'stim'
        try
            ST = loadcb(cellid,'STIMSPIKES');
            TE = loadcb(cellid,'StimEvents');
        catch
            out = [NaN NaN NaN NaN];
        end
    case 'behav'
        try
           ST = loadcb(cellid,'EVENTSPIKES');
            TE = loadcb(cellid,'TrialEvents');
        catch
            out = [NaN NaN NaN NaN];
        end
end
epoch_pos = findcellstr(ST.epochs(:,1),EpochName);

if (epoch_pos == 0)
    error('Epoch name not found');
end

RATE = ST.epoch_rates{epoch_pos};

try
    ValidTrials = selecttrial(TE,'OdorPokeValid & WaterPokeValid ');
catch
    ValidTrials=1:length(TE.TrialStart);
end

if lookback == 1,
    posOMIT      = intersect(find(TE.PreviousRewardOmit==1),ValidTrials);
    posREWARD       = intersect(find(TE.PreviousRewardOmit==0),ValidTrials);
else
    posOMIT      = intersect(find(TE.RewardOmit==1),ValidTrials);
    posREWARD       = intersect(find(TE.RewardOmit==0),ValidTrials);
end
[Domi,  Pomi] = rocarea(RATE(posREWARD),RATE(posOMIT),'boot',nboot,'scale');

% Here is how you'd do side selectivity:
%
%posLEFT       = intersect(find(TE.ChoiceLeft),ValidTrials);
%posRIGHT      = intersect(find(TE.ChoiceRight),ValidTrials);
%[Dside, Pside]   =  rocarea(RATE(posLEFT),RATE(posRIGHT),'boot',nboot,'scale');

out = [Domi Pomi length(posOMIT) length(posREWARD)];
