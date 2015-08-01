function out = calc_maxrate_perf0(cellid,EpochName);
% 
% CALC_MAXRATE_PREF
%
% How much could the rat improve by switching choices above a certain rate
% threshold?
%
%
% [ImprovRatio, NumTrials2Switch, AccuBelow] = calc_choiceswitch_improv(cellid,EpochName)
% 

%AK 3/07

% if strcmpi(cellid,'default')
%     out = [strcat({'Dout_','Pout_'},{EpochName,EpochName})];
%     return;
% end

TE = loadcb(cellid,'Events');
ST = loadcb(cellid,'EVENTSPIKES');

epoch_pos = findcellstr(ST.epochs(:,1),EpochName);

if (epoch_pos == 0)
    error('Epoch name not found');
end

RATE = ST.epoch_rates{epoch_pos};
ValidTrials = selecttrial(TE,'OdorPokeValid & WaterPokeValid & OdorPairID == 1');
%
RATE = RATE(ValidTrials);
CORRECT = TE.Correct(ValidTrials);
%

C90=CORRECT(find(RATE>prctile(RATE,90)));
accu90 = nanmean([C90 NaN]);

C95=CORRECT(find(RATE>prctile(RATE,95)));
accu95 = nanmean([C95 NaN]);

[m,v] =binostat(1, [accu90 accu95]);
STE = sqrt(v)./ max(1,[sqrt(length(C90)-1) sqrt(length(C90)-1)]);

%out = [accu90 accu95 STE];

out = STE;