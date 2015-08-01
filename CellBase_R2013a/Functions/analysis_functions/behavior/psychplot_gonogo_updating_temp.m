function psychplot_gonogo_updating(cbID,varargin)
%PSYCHPLOT_PAVLOVIAN2   Plot performance.
%   PSYCHPLOT_PAVLOVIAN2(CELLID) calculates performance and reaction time
%   variables and plots lick rasters. The session corresponding to the cell
%   given in CELLID is selected.
%
%   PSYCHPLOT_PAVLOVIAN2(SESSIONID) calculates performance and reaction
%   time variables for the session passed in SESSIONID. SESSIONID should be
%   a 1-by-2 cell array with the animal ID and the session ID.
%
%   Additional input parameter-value pairs, with default values:
%       'valid_trials', 'all' - restrict the analysis to a set of trials
%       'lastresponsestop', false - exclude trials after the last response
%           of the animal
%       'display', true - controls plotting
%
%   Output:
%       PerformanceR1 - lick probability for water in high prob. reward trials
%       PerformanceP1 - lick probability for airpuff in high prob. reward trials
%       PerformanceO1 - lick probability for omissions in high prob. reward trials
%       PerformanceR2 - lick probability for water in low prob. reward trials
%       PerformanceP2 - lick probability for airpuff in low prob. reward trials
%       PerformanceO2 - lick probability for omissions in low prob. reward trials
%       AnticipationR1 - lick probability in delay before water in high prob. reward trials
%       AnticipationP1 - lick probability in delay before airpuff in high prob. reward trials
%       AnticipationO1 - lick probability in delay before omissions in high prob. reward trials
%       AnticipationR2 - lick probability in delay before airpuff in low prob. reward trials
%       AnticipationP2 - lick probability in delay before airpuff in low prob. reward trials
%       AnticipationO2 - lick probability in delay before omissions in low prob. reward trials
%       ReactionTimeR1 - reaction time to water in high prob. reward trials
%       ReactionTimeP1 - reaction time to airpuff in high prob. reward trials
%       ReactionTimeO1 - reaction time to omissions in high prob. reward trials
%       ReactionTimeR2 - reaction time to water in low prob. reward trials
%       ReactionTimeP2 - reaction time to airpuff in low prob. reward trials
%       ReactionTimeO2 - reaction time to omissions in low prob. reward trials
%
%   See also PSYCHPLOT_PAVLOVIAN.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   13-May-2014

%   Edit log: BH 5/13/14

% Input arguments
prs = inputParser;
addRequired(prs,'cbID',@(s)iscellid(s)|issessionid(s))   % cell ID or session ID
% addRequired(prs,'cbID')   % cell ID or session ID
addParamValue(prs,'valid_trials','all',...
    @(s)(ischar(s)&isequal(s,'all'))|isnumeric(s))   % valid trials
addParamValue(prs,'lastresponsestop',false,@(s)islogical(s)|ismember(s,[0 1]))   % exclude trials after the last response
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
hit_inx = find(TE.Hit==1);
fa_inx = find(TE.FalseAlarm==1);
[GPhp NGPhp GORThp NGRThp] = psychplot_gonogo({animalID sessionID},'display',false,...
    'valid_trials',intersect(g.valid_trials,hit_inx+1));
[GPhm NGPhm GORThm NGRThm] = psychplot_gonogo({animalID sessionID},'display',false,...
    'valid_trials',intersect(g.valid_trials,hit_inx-1));
[GPfp NGPfp GORTfp NGRTfp] = psychplot_gonogo({animalID sessionID},'display',false,...
    'valid_trials',intersect(g.valid_trials,fa_inx+1));
[GPfm NGPfm GORTfm NGRTfm] = psychplot_gonogo({animalID sessionID},'display',false,...
    'valid_trials',intersect(g.valid_trials,fa_inx-1));

% Plot
if g.display
    
    figure
    viewlick({animalID sessionID},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#FeedbackID','window',[-5 5])
    maximize_figure(gcf)
    
    figure
    viewlick({animalID sessionID},'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#StimulusID','window',[-5 5])
    maximize_figure(gcf)
    
    figure
    viewlick({animalID sessionID},'TriggerName','DeliverFeedback','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#FeedbackID','window',[-5 5])
    maximize_figure(gcf)
    
    rasterpsth({animalID sessionID})
    
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
    bar(2,PerformanceP1,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(3,PerformanceO1,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(4.5,PerformanceR2,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    bar(5.5,PerformanceP2,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(6.5,PerformanceO2,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    set(gca,'XTick',[1 2.5 4 5.5],'XTickLabel',{'R1' 'R2' 'P1' 'P2'})
    box off
    title('Performance')
    
    % Anticipation plot
    subplot(232)
    bar(1,AnticipationR1,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    hold on
    bar(2,AnticipationP1,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(3,AnticipationO1,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(4.5,AnticipationR2,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    bar(5.5,AnticipationP2,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(6.5,AnticipationO2,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    set(gca,'XTick',[1 2.5 4 5.5],'XTickLabel',{'R1' 'R2' 'P1' 'P2'})
    box off
    title('Anticipation')
    
    % ReactionTime plot
    subplot(233)
    bar(1,ReactionTimeR1,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    hold on
    bar(2,ReactionTimeP1,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(3,ReactionTimeO1,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
    bar(4.5,ReactionTimeR2,'BarWidth',0.5,'EdgeColor',rwcolor,'FaceColor','w','LineWidth',2)
    bar(5.5,ReactionTimeP2,'BarWidth',0.5,'EdgeColor',pncolor,'FaceColor','w','LineWidth',2)
    bar(6.5,ReactionTimeO2,'BarWidth',0.5,'EdgeColor','k','FaceColor','w','LineWidth',2)
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

% -------------------------------------------------------------------------
function rasterpsth(cellid)

% Time window
wn = [-5 5];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector

% Calcualte PSTH
[psth, spsth, spsth_se, ~, spt] = ...
    ultimate_psth(cellid,'lick','StimulusOn',wn,...
    'dt',dt,'sigma',0.02,'parts','#FeedbackID','isadaptive',2,...
    'maxtrialno',Inf,'first_event','FirstITIBegins','last_event','TrialEnd');

% Lick raster
H = figure;
viewlick(cellid,'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav',...
    'ShowEvents',{{'FirstITIBegins' 'StimulusOn' 'StimulusOff' 'DelayEnd' 'TrialEnd'}},...
    'Partitions','#FeedbackID','window',[-5 5])
maximize_figure(H)

% Replace PSTH
L = findobj(allchild(gca),'Type','line');   % lines in the legend
clr = flipud(get(L,'Color'));
cla
for k = 1:size(psth,1)
    errorshade(time,spsth(k,:),spsth_se(k,:),'LineWidth',2,...
        'LineColor',clr{k},'ShadeColor',clr{k})
    hold on
end
axis tight
set(gcf,'Renderer','OpenGL')

figure
errorshade(time,spsth(1,:),spsth_se(1,:),'LineWidth',2,...
        'LineColor',clr{1},'ShadeColor',clr{1})

% Save
% fnm = [resdir cellidt '_' alignevent '_rasterPSTH.jpg'];   % save
% saveas(H,fnm)
% if issave
%     set(H,'PaperPositionMode','auto')
%     set(H,'InvertHardcopy','off')
%     print(H,'-djpeg',fnm)
% end