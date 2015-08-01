function out = calc_choiceswitch_improvement(cellid,EpochName);
% 
% CALC_CHOICESWITCH_IMPROV
%
% How much could the rat improve by switching choices above a certain rate
% threshold?
%
%
% [ImprovRatio, NumTrials2Switch, AccuBelow] = calc_choiceswitch_improv(cellid,EpochName)
% 

%AK 3/07

if strcmpi(cellid,'default')
    out = [strcat({'Dout_','Pout_'},{EpochName,EpochName})];
    return;
end

TE = loadcb(cellid,'Events');
ST = loadcb(cellid,'EVENTSPIKES');

epoch_pos = findcellstr(ST.epochs(:,1),EpochName);

if (epoch_pos == 0)
    error('Epoch name not found');
end

RATE = ST.epoch_rates{epoch_pos};
try
    ValidTrials = selecttrial(TE,'OdorPokeValid & WaterPokeValid & OdorPairID == 1');
catch
    ValidTrials=ones(1,length(TE));
end
%
RATE = RATE(ValidTrials);
CORRECT = TE.Correct(ValidTrials);
%
URATE = sort(unique(RATE)); URATE=URATE(end:-1:1);

ACCU(1)=NaN;
STE(1)=NaN;
for iR=2:length(URATE)
    posRATE = find(RATE>URATE(iR));
    ACCU(iR) = mean(TE.Correct(posRATE));
    [m,v] =binostat(1, ACCU(iR)); 
    STE(iR) = sqrt(v)./max(1,sqrt(length(posRATE)-1));
end

posIMPROV=find(ACCU+STE<0.5);                   %rate threshold for improvement;
if ~isempty(posIMPROV)
    posTRIALS = find(RATE>min(URATE(posIMPROV)));
    AccuBelow=mean(CORRECT(posTRIALS));
    NumTrials2Switch = length(find(CORRECT(posTRIALS)==0)) - length(find(CORRECT(posTRIALS)==1));
else
    NumTrials2Switch = 0;
    AccuBelow = NaN;
end
TotalTrials = length(ValidTrials);
TotalCorrect = nansum(CORRECT);

ImprovRatio=(TotalCorrect+NumTrials2Switch)*TotalTrials/(TotalTrials*TotalCorrect);

out = [ImprovRatio, NumTrials2Switch, AccuBelow];

