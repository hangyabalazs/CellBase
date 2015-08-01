function out = calc_selectivity_side(cellid,EpochName,nboot,eventtype)
% 
% CALC_SELECTIVITY_SIDE
%
% Calculates the area under the ROC measure for side selectivity in the
% requested epoch and a bootsrapped P-value.
%
% calc_selectivity_side(cellid,EpochName,nboot)

% AK 10/06 SPR 2010-08-04

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

if strcmpi(cellid,'default')
    out = [strcat({'Dside_','Pside_'},{EpochName,EpochName})];
    return;
end

epoch_pos = findcellstr(ST.epochs(:,1),EpochName);

if (epoch_pos == 0)
    error('Epoch name not found');
end

RATE = ST.epoch_rates{epoch_pos};

try
    ValidTrials = selecttrial(TE,'OdorPokeValid & WaterPokeValid ');
    posLEFT        = intersect(find(TE.ChoiceLeft),ValidTrials);
    posRIGHT       = intersect(find(TE.ChoiceRight),ValidTrials);
catch
    ValidTrials=1:length(TE.TrialStart);
    posLEFT = intersect(find(TE.SoundID==1),ValidTrials);
    posRIGHT = intersect(find(TE.SoundID==2),ValidTrials);
end

[Dside, Pside] = rocarea(RATE(posLEFT),RATE(posRIGHT),'boot',nboot,'scale');

out = [Dside Pside];
