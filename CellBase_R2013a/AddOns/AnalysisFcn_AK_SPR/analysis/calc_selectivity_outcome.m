function out = calc_selectivity_outcome(cellid,EpochName,nboot,eventtype);
% 
% CALC_SELECTIVITY_OUTCOME
%
% Calculates the area under the ROC measure for outcome selectivity in the
% requested epoch and a bootsrapped P-value.
%
% calc_selectivity_outcome(cellid,EpochName,nboot)

%AK 10/06

if strcmpi(cellid,'default')
    out = [strcat({'Dout_','Pout_'},{EpochName,EpochName})];
    return;
end

if nargin<4,
    eventtype='behav';
end
switch eventtype
    case 'stim'
        ST = loadcb(cellid,'STIMSPIKES');
        TE = loadcb(cellid,'StimEvents');
    case 'behav'
        ST = loadcb(cellid,'EVENTSPIKES');
        TE = loadcb(cellid,'TrialEvents');
end
epoch_pos = findcellstr(ST.epochs(:,1),EpochName);

if (epoch_pos == 0)
    error('Epoch name not found');
end

RATE = ST.epoch_rates{epoch_pos};

try
    ValidTrials = selecttrial(TE,'OdorPokeValid & WaterPokeValid ');
catch
    ValidTrials=ones(1,length(TE));
end

posERROR      = intersect(find(TE.Error),ValidTrials);
posCORR       = intersect(find(TE.Correct),ValidTrials);
[Dout,  Pout] = rocarea(RATE(posCORR),RATE(posERROR),'boot',nboot,'scale');

% Here is how you'd do side selectivity:
%
%posLEFT       = intersect(find(TE.ChoiceLeft),ValidTrials);
%posRIGHT      = intersect(find(TE.ChoiceRight),ValidTrials);
%[Dside, Pside]   =  rocarea(RATE(posLEFT),RATE(posRIGHT),'boot',nboot,'scale');

out = [Dout Pout];
