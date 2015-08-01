function [P, PM, iG] = get_performance(te,varargin)
%
%   [OdorRatio ProbLeft SD  SE NumStim RewProb]
%

STIM = unique(te.OdorRatio);
NumOdorPairs = length(unique(te.OdorPairID));
PM = zeros(NumOdorPairs,length(STIM),6)*NaN;
PM(:,:,1)=repmat(STIM,NumOdorPairs,1);

% N=0;
% for iG=1:NumOdorPairs; 
%      N = N+length(unique(te.OdorConc(find(te.OdorPairID==iG))));
% end
if nargin > 1
    pos2USE = varargin{1};
else
    pos2USE = 1:length(te.OdorRatio);
end

for iG=1:NumOdorPairs;
    posG = intersect(find(te.OdorPairID==iG),pos2USE);
    stim = unique(te.OdorRatio(posG));
    clear psycho;
    for iS = 1:length(stim)
        posS = intersect(find(te.OdorRatio==stim(iS)),posG);
        stimNUM = length(posS);
        psycho(iS,1) = stim(iS);
        psycho(iS,2) = mean(te.ChoiceLeft(posS));
        psycho(iS,3) = std(te.ChoiceLeft(posS));
        psycho(iS,4) = psycho(iS,3)/sqrt(stimNUM-1);
        psycho(iS,5) = stimNUM;
        psycho(iS,6) = nanmean(te.RewardProb(posS));
        iSTIM = find(STIM==stim(iS));
       % PM(iG,iSTIM,1) = stim(iS);
        PM(iG,iSTIM,2) = mean(te.ChoiceLeft(posS));
        PM(iG,iSTIM,3) = std(te.ChoiceLeft(posS));
        PM(iG,iSTIM,4) = psycho(iS,3)/sqrt(stimNUM-1);
        PM(iG,iSTIM,5) = stimNUM;
        PM(iG,iSTIM,6) = nanmean(te.RewardProb(posS));
    end %iS
    P(iG).psycho = psycho;
end %iG


