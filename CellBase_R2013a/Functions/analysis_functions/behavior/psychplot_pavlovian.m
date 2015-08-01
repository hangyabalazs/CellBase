function [PerformanceR1 PerformanceR1O PerformanceR2 PerformanceR2O ...
    PerformanceP1 PerformanceP1O PerformanceP2 PerformanceP2O ...
    AnticipationR1 AnticipationR1O AnticipationR2 AnticipationR2O ...
    AnticipationP1 AnticipationP1O AnticipationP2 AnticipationP2O ...
    ReactionTimeR1 ReactionTimeR1O ReactionTimeR2 ReactionTimeR2O ...
    ReactionTimeP1 ReactionTimeP1O ReactionTimeP2 ReactionTimeP2O] = ...
    psychplot_pavlovian(cbID,varargin)
%PSYCHPLOT_PAVLOVIAN   Plot performance.
%   PSYCHPLOT_PAVLOVIAN(CELLID) calculates performance and reaction time
%   variables and plots lick rasters. The session corresponding to the cell
%   given in CELLID is selected.
%
%   PSYCHPLOT_PAVLOVIAN(SESSIONID) calculates performance and reaction time
%   variables for the session passed in SESSIONID. SESSIONID should be a
%   1-by-2 cell array with the animal ID and the session ID.
%
%   Additional input parameter-value pairs, with default values:
%       'valid_trials', 'all' - restrict the analysis to a set of trials
%       'lastresponsestop', false - exclude trials after the last response
%           of the animal
%       'display', true - controls plotting
%
%   Output:
%       PerformanceR1 - lick probability for water in high prob. reward trials
%       PerformanceR1O - lick probability for omissions in high prob. reward trials
%       PerformanceR2 - lick probability for water in low prob. reward trials
%       PerformanceR2O - lick probability for omissions in low prob. reward trials
%       PerformanceP1 - lick probability for airpuff in high prob. punishment trials
%       PerformanceP1O - lick probability for omissions in high prob. punishment trials
%       PerformanceP2 - lick probability for airpuff in low prob. punishment trials
%       PerformanceP2O - lick probability for omissions in low prob. punishment trials
%       AnticipationR1 - lick probability in delay before water in high prob. reward trials
%       AnticipationR1O - lick probability in delay before omissions in high prob. reward trials
%       AnticipationR2 - lick probability in delay before water in low prob. reward trials
%       AnticipationR2O - lick probability in delay before omissions in low prob. reward trials
%       AnticipationP1 - lick probability in delay before airpuff in high prob. punishment trials
%       AnticipationP1O - lick probability in delay before omissions in high prob. punishment trials
%       AnticipationP2 - lick probability in delay before airpuff in low prob. punishment trials
%       AnticipationP2O - lick probability in delay before omissions in low prob. punishment trials
%       ReactionTimeR1 - reaction time to water in high prob. reward trials
%       ReactionTimeR1O - reaction time to omissions in high prob. reward trials
%       ReactionTimeR2 - reaction time to water in low prob. reward trials
%       ReactionTimeR2O - reaction time to omissions in low prob. reward trials
%       ReactionTimeP1 - reaction time to airpuff in high prob. punishment trials
%       ReactionTimeP1O - reaction time to omissions in high prob. punishment trials
%       ReactionTimeP2 - reaction time to airpuff in low prob. punishment trials
%       ReactionTimeP2O - reaction time to omissions in low prob. punishment trials
%
%   See also PSYCHPLOT_GONOGO.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   1-May-2014

%   Edit log: BH 5/1/14

% Input arguments
prs = inputParser;
% addRequired(prs,'cbID',@(s)iscellid(s)|issessionid(s))   % cell ID or session ID
addRequired(prs,'cbID')   % cell ID or session ID
addParamValue(prs,'valid_trials','all',...
    @(s)(ischar(s)&isequal(s,'all'))|isnumeric(s))   % valid trials
addParamValue(prs,'lastresponsestop',false,@(s)islogical(s)|ismember(s,[0 1]))   % control ploting
addParamValue(prs,'display',true,@(s)islogical(s)|ismember(s,[0 1]))   % control ploting
parse(prs,cbID,varargin{:})
g = prs.Results;
if (ischar(g.cbID) || (iscell(cbID) && length(g.cbID)==1)) && iscellid(g.cbID)   % cell ID was passed
    [animalID sessionID] = cellid2tags(g.cbID);
else   % session ID was passed
    animalID = g.cbID{1,1};
    sessionID = g.cbID{1,2};
end

% Load trial events (unsynchronized)
cbdir = getpref('cellbase','datapath');
datapath = fullfile(cbdir,animalID,sessionID,'TE.mat');
TE = load(datapath);
if isequal(g.valid_trials,'all')   % valid trials
    g.valid_trials = 1:length(TE.TrialStart)-1;  % last trial may be incomplete
end
if g.lastresponsestop   % drop trials after last hit
    lastresp = find(TE.Hit==1|TE.FalseAlarm==1,1,'last');
    g.valid_trials(g.valid_trials>lastresp) = [];
end
TE = filterTE(TE,g.valid_trials);

% Performance
PerformanceR1 = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.R1Trial)) / ...
    sum(nan2zero(TE.R1Trial));
PerformanceR2 = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.R2Trial)) / ...
    sum(nan2zero(TE.R2Trial));
PerformanceR1O = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.R1OTrial)) / ...
    sum(nan2zero(TE.R1OTrial));
PerformanceR2O = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.R2OTrial)) / ...
    sum(nan2zero(TE.R2OTrial));

PerformanceP1 = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.P1Trial)) / ...
    sum(nan2zero(TE.P1Trial));
PerformanceP2 = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.P2Trial)) / ...
    sum(nan2zero(TE.P2Trial));
PerformanceP1O = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.P1OTrial)) / ...
    sum(nan2zero(TE.P1OTrial));
PerformanceP2O = sum(nan2zero(TE.FeedbackLick)&nan2zero(TE.P2OTrial)) / ...
    sum(nan2zero(TE.P2OTrial));

% Delay licks
AnticipationR1 = sum(nan2zero(TE.DelayLick)&nan2zero(TE.R1Trial)) / ...
    sum(nan2zero(TE.R1Trial));
AnticipationR2 = sum(nan2zero(TE.DelayLick)&nan2zero(TE.R2Trial)) / ...
    sum(nan2zero(TE.R2Trial));
AnticipationR1O = sum(nan2zero(TE.DelayLick)&nan2zero(TE.R1OTrial)) / ...
    sum(nan2zero(TE.R1OTrial));
AnticipationR2O = sum(nan2zero(TE.DelayLick)&nan2zero(TE.R2OTrial)) / ...
    sum(nan2zero(TE.R2OTrial));

AnticipationP1 = sum(nan2zero(TE.DelayLick)&nan2zero(TE.P1Trial)) / ...
    sum(nan2zero(TE.P1Trial));
AnticipationP2 = sum(nan2zero(TE.DelayLick)&nan2zero(TE.P2Trial)) / ...
    sum(nan2zero(TE.P2Trial));
AnticipationP1O = sum(nan2zero(TE.DelayLick)&nan2zero(TE.P1OTrial)) / ...
    sum(nan2zero(TE.P1OTrial));
AnticipationP2O = sum(nan2zero(TE.DelayLick)&nan2zero(TE.P2OTrial)) / ...
    sum(nan2zero(TE.P2OTrial));

% Reaction time
ReactionTimeR1 = nanmedian(TE.Reward1RT);
ReactionTimeR2 = nanmedian(TE.Reward2RT);
ReactionTimeR1O = nanmedian(TE.Reward1OmissionRT);
ReactionTimeR2O = nanmedian(TE.Reward2OmissionRT);
ReactionTimeP1 = nanmedian(TE.Punish1RT);
ReactionTimeP2 = nanmedian(TE.Punish2RT);
ReactionTimeP1O = nanmedian(TE.Punish1OmissionRT);
ReactionTimeP2O = nanmedian(TE.Punish2OmissionRT);

% Plot
if g.display
    
    figure
    viewlick({animalID sessionID},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#FeedbackID','window',[-5 5])
    maximize_figure(gcf)
    
    figure
    viewlick({animalID sessionID},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#StimulusID','window',[-5 5])
    maximize_figure(gcf)
    
    % Plot parameters
    rwcolor = [0 0.8 0];   % color for reward performace
    pncolor = [0.8 0 0];   % color for punishment performance
    rtcolor = 'k';   % color for reaction time
    ms = 8;   % marker size
    lw = 2;   % line width
    
    % Open figure
    H = figure;
    set(gcf,'DefaultAxesFontName','Arial')
    set(gcf,'DefaultAxesTickDir','out')
    set(gcf,'DefaultAxesBox','off')
    
    % Performance plot
    subplot(231)
    bar(1,PerformanceR1,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    hold on
    bar(1.5,PerformanceR1O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(2.5,PerformanceR2,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    bar(3,PerformanceR2O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(4,PerformanceP1,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(4.5,PerformanceP1O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(5.5,PerformanceP2,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(6,PerformanceP2O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    set(gca,'XTick',[1 2.5 4 5.5],'XTickLabel',{'R1' 'R2' 'P1' 'P2'})
    box off
    title('Performance')
    
    % Anticipation plot
    subplot(232)
    bar(1,AnticipationR1,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    hold on
    bar(1.5,AnticipationR1O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(2.5,AnticipationR2,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    bar(3,AnticipationR2O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(4,AnticipationP1,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(4.5,AnticipationP1O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(5.5,AnticipationP2,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(6,AnticipationP2O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    set(gca,'XTick',[1 2.5 4 5.5],'XTickLabel',{'R1' 'R2' 'P1' 'P2'})
    box off
    title('Anticipation')
    
    % ReactionTime plot
    subplot(233)
    bar(1,ReactionTimeR1,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    hold on
    bar(1.5,ReactionTimeR1O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(2.5,ReactionTimeR2,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    bar(3,ReactionTimeR2O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(4,ReactionTimeP1,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(4.5,ReactionTimeP1O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(5.5,ReactionTimeP2,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(6,ReactionTimeP2O,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    set(gca,'XTick',[1 2.5 4 5.5],'XTickLabel',{'R1' 'R2' 'P1' 'P2'})
    box off
    title('Reaction Time')
    maximize_figure(H)
    
    % Sort the trials into non-overlapping blocks
%     [Bl Ia Ib] = unique(TE.BlockNum(1:end-1));
%     [bGoPerformance bNoGoPerformance] = deal(nan(1,length(Bl)));  % performance per block
%     for ibl = 1:length(Bl),
%         bGoPerformance (ibl) = nansum(TE.Hit(Ib==ibl))/(nansum(TE.Hit(Ib==ibl))+nansum(TE.Miss(Ib==ibl)));
%         bNoGoPerformance(ibl) = nansum(TE.FalseAlarm(Ib==ibl))/(nansum(TE.FalseAlarm(Ib==ibl))+nansum(TE.CorrectRejection(Ib==ibl)));
%     end
%     
%     % Plot block-wise performance
%     subplot(2,3,4:6)
%     plot(bGoPerformance,'o-','Color',gocolor);
%     hold on
%     plot(bNoGoPerformance,'o-','Color',nogocolor);
%     plot((bGoPerformance - bNoGoPerformance),'ks-')
%     legend({'%Hit','%FA','d'''},'Location','best')
end

% Supress output
if nargout < 1
    clear
end

% -------------------------------------------------------------------------
function TE = filterTE(TE,valid_trials)

fnm = fieldnames(TE);
for k = 1:length(fnm)
    TE.(fnm{k}) = TE.(fnm{k})(valid_trials);
end