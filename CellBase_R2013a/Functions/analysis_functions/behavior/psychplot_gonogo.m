function [GoPerformance NoGoPerformance GoReactionTime NoGoReactionTime] = ...
    psychplot_gonogo(cbID,varargin)
%PSYCHPLOT_GONOGO   Plot psychometric performance.
%   [GOPERF NGPERF GORT NGRT] = PSYCHPLOT_GONOGO(CELLID) calculates go and 
%   no-go performance (GOPERF and NGPERF) as well as go and no-go reaction
%   time (GORT and NGRT) for all different stimulus intensities in the
%   auditory go/no-go task. The session corresponding to the cell given in
%   CELLID is selected.
%
%   [GOPERF NGPERF GORT NGRT] = PSYCHPLOT_GONOGO(SESSIONID) calculates
%   performance and reaction time for the session passed in SESSIONID.
%   SESSIONID should be a 1-by-2 cell array with the animal ID and the
%   session ID.
%
%   Additional input parameter-value pairs, with default values:
%       'valid_trials', 'all' - restrict the analysis to a set of trials
%       'lastresponsestop', false - exclude trials after the last response
%           of the animal
%       'display', true - controls plotting; performance, median reaction
%           time, reaction time cumulative density functions, difference
%           between go and no-go performace and temporal evolution of
%           blockwise performance are plotted
%
%   See also PSYCH_GONOGO.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Sept-2013

%   Edit log: BH 9/27/13

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
if (ischar(g.cbID) || (iscell(cbID) && length(g.cbID)==1)) && iscellid(g.cbID)   % cell ID was passed   % cell ID was passed
    [animalID sessionID] = cellid2tags(g.cbID);
else   % session ID was passed
    animalID = g.cbID{1,1};
    sessionID = g.cbID{1,2};
end

% Load trial events (unsynchronized)
cbdir = getpref('cellbase','datapath');
datapath = fullfile(cbdir,animalID,sessionID,'TE.mat');
TE = load(datapath);
if isequal(whichcb,'human_susattn')   % different sound coding on the human rig
    TE.StimulusDuration = log(TE.SoundIntensity);
end
if isequal(g.valid_trials,'all')   % valid trials
    g.valid_trials = 1:length(TE.TrialStart)-1;  % last trial may be incomplete
end
if g.lastresponsestop   % drop trials after last hit
    lastresp = find(TE.Hit==1|TE.FalseAlarm==1,1,'last');
    g.valid_trials(g.valid_trials>lastresp) = [];
end
TE = filterTE(TE,g.valid_trials);

% Performance for all sound intensities
[SI Ia Ib] = unique(TE.StimulusDuration(1:end-1));   % find unique sound intensities
NumInt = length(Ia);   % number of different stimuli
[GoPerformance NoGoPerformance GoReactionTime NoGoReactionTime] = deal(nan(1,NumInt));
for iS = 1:NumInt
    
    % Go performance
    GoPerformance(iS) = nansum(TE.Hit(Ib==iS)) / ...
        (nansum(TE.Hit(Ib==iS)) + nansum(TE.Miss(Ib==iS)));
    GoReactionTime(iS) = nanmedian(TE.GoRT(Ib==iS));
    
    % No-go performance
    NoGoPerformance(iS) = nansum(TE.FalseAlarm(Ib==iS)) / ...
        (nansum(TE.CorrectRejection(Ib==iS)) + nansum(TE.FalseAlarm(Ib==iS)));
    NoGoReactionTime(iS) = nanmedian(TE.NoGoRT(Ib==iS));
end

% Plot
if g.display
    
    % Plot parameters
    gocolor = [0 0.8 0];   % color for go performace
    nogocolor = [0.8 0 0];   % color for no-go performance
    rtcolor = 'k';   % color for reaction time
    dcolor = 'k';   % solor for pseudo-d-prime
    ms = 8;   % marker size
    lw = 2;   % line width
    
    % Open figure
    H = figure;
    set(gcf,'DefaultAxesFontName','Arial')
    set(gcf,'DefaultAxesTickDir','out')
    set(gcf,'DefaultAxesBox','off')
    
    % Performance plot
    subplot(251)
    plot(SI,GoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
    hold on
    plot(SI,NoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
    ylim([0 1])
    box off
    title('Performance')
    
    subplot(252)
    plot(SI,GoReactionTime,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',rtcolor,'Color',rtcolor,'MarkerSize',ms,'LineWidth',lw);
    hold on
    plot(SI,NoGoReactionTime,'LineStyle','--','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',rtcolor,'Color',rtcolor,'MarkerSize',ms,'LineWidth',lw);
    box off
    title('Reaction time')
    xlabel('Sound Pressure Level (dB)')
    
    % Reaction time distributions
    bno = 100;  % bin number
    edges = linspace(nanmin(TE.ReactionTime),nanmax(TE.ReactionTime),bno);
    shadeness = linspace(0.9,0.1,NumInt);
    subplot(253)
    hold on
    for iSub = 1:NumInt   % Go reaction time CDF
        N = histc(TE.GoRT(TE.StimulusDuration==SI(iSub)),edges);
        N = [0 N(1:end-1)];
        stairs(edges,cumsum(N)./max(cumsum(N)),'Color',[0 shadeness(iSub) 0],'LineWidth',lw);
        axis tight
    end
    subplot(254)
    hold on
    for iSub = 1:NumInt   % No-go reaction time CDF
        N = histc(TE.NoGoRT(TE.StimulusDuration==SI(iSub)),edges);
        N = [0 N(1:end-1)];
        stairs(edges,cumsum(N)./max(cumsum(N)),'Color',[shadeness(iSub) 0 0],'LineWidth',lw);
        axis tight
    end
    
    % Pseudo-d-prime
    subplot(255)
    hold on
    plot(SI,GoPerformance-NoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',dcolor,'Color',dcolor,'MarkerSize',ms,'LineWidth',lw);
    ylim([-1 1])
    title('d''')
    fstamp([animalID ' ' sessionID])
    
    % Sort the trials into non-overlapping blocks
    [Bl Ia Ib] = unique(TE.BlockNum(1:end-1));
    [bGoPerformance bNoGoPerformance] = deal(nan(1,length(Bl)));  % performance per block
    for ibl = 1:length(Bl),
        bGoPerformance (ibl) = nansum(TE.Hit(Ib==ibl))/(nansum(TE.Hit(Ib==ibl))+nansum(TE.Miss(Ib==ibl)));
        bNoGoPerformance(ibl) = nansum(TE.FalseAlarm(Ib==ibl))/(nansum(TE.FalseAlarm(Ib==ibl))+nansum(TE.CorrectRejection(Ib==ibl)));
    end
    
    % Plot block-wise performance
    subplot(2,3,4:6)
    plot(bGoPerformance,'o-','Color',gocolor);
    hold on
    plot(bNoGoPerformance,'o-','Color',nogocolor);
    plot((bGoPerformance - bNoGoPerformance),'ks-')
    legend({'%Hit','%FA','d'''},'Location','best')
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