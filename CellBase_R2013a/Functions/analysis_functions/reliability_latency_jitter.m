function [R L J B M lim1 lim2 spno H] = reliability_latency_jitter(cellid,varargin)
%RELIABILITY_LATENCY_JITTER   Reliability, latency and jitter of spiking.
%   [R L J] = RELIABILITY_LATENCY_JITTER(CELLID) calculates the latency of
%   spiking of a cell (CELLID) after an event (default, 'PulseOn') (L). It
%   also returns the standard deviation of spike times after finding the
%   period of evoked spikes using ULTIMATE_PSTH (J). Reliability of evoking
%   spikes (proportion of events followed by evoked spikes) is returned in
%   R.
%
%   [R L J B M] = RELIABILITY_LATENCY_JITTER(CELLID,...) also returns
%   baseline firing rate (B) and maximal firing rate in the test window
%   (M).
%
%   [R L J B M A1 A2] = RELIABILITY_LATENCY_JITTER(CELLID,...) returns the
%   start and end point of detected stimulation period (see
%   ULTIMATE_PSTH and FINDSTIMPERIOD).
%
%   [R L J B M A1 A2 SPNO] = RELIABILITY_LATENCY_JITTER(CELLID,...) returns
%   the number of spikes detected for each valid trial (see EFFICIENCY).
%
%   [R L J B M A1 A2 SPNO H] = RELIABILITY_LATENCY_JITTER(CELLID,...)
%   returns figure handles for PSTH and raster plot in the structure H (if
%   'display' is set to true).
%
%   Optional parameter-value pairs, with default values:
%       'event_type', 'stim' - type of the aligning event, 'stim' (stimulus
%           events) or 'trial' (trial events)
%       'event', 'PulseOn' - aligning event
%       'event_filter', 'none' - filter trials; see FILTERTRIALS for 
%           implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%       'window', [-0.005 0.01] - time window relative to the event, in
%           seconds; for determining the range of the PSTH, see
%           ULTIMATE_PSTH
%       'isadaptive', 1 - 0, classic PSTH algorithm is applied; 1, adaptive
%           PSTH is calculated (see APSTH); 2, 'doubly adaptive' PSTH
%           algorithm is used (see DAPSTH)
%   	'baselinewin', [-0.005 0] - limits of the baseline window for PSTH 
%           statistics (see PSTH_STATS), time relative to 0 in seconds
%   	'testwin', [0 0.01] - limits of the test window for PSTH statistics
%           (see PSTH_STATS), time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals in PSTH_STATS; in
%           proportion of the peak-baseline difference (see PSTH_STATS); it
%           determines the window within which spikes will be considered
%           'evoked'
%       'jitterdefinition', 'all' - control how jitter is defined: 
%           'all' - SD of spike times of all spikes
%           'burst' - only first spikes are included in each trial
%       'display', false - controls plotting
%
%   Examples:
%   [E L J] = reliability_latency_jitter(cellid,'event','BurstOn');
%
%   fi = struct('BurstNPulse',20);
%   [E I J] = reliability_latency_jitter(cellid,...
%       'event_filter','BurstNPulse_maxPower','filterinput',fi);
%
%   [E_burston L_burston J_burston B_burston M_burston Astart Aend] = ...
%        reliability_latency_jitter(cellid,'event','PulseOn');
%
%   See also ULTIMATE_PSTH, EXTRACTSEGSPIKES and EFFICIENCY.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   15-June-2013

%   Edit log: BH 7/15/13

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'event_type','stim',...
    @(s)ischar(s)&&(strncmp(s,'stim',4)||strncmp(s,'trial',4)))   % event type - stim or trial
addParamValue(prs,'event','PulseOn',@ischar)   % reference event
addParamValue(prs,'event_filter','none',@ischar)   % event filter
addParamValue(prs,'filterinput',[])   % some filters need additional input
addParamValue(prs,'window',[-0.005 0.01],...
    @(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'isadaptive',1,@(s)islogical(s)|ismember(s,[0 1 2]))   % use adaptive PSTH algorithm
addParamValue(prs,'baselinewin',[-0.005 0],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParamValue(prs,'testwin',[0 0.01],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParamValue(prs,'relative_threshold',0.5,...
    @(s)isnumeric(s)&isequal(length(s),1)&s>=-1&s<=1)  % relative threshold for peak detection, see ULTIMATE_PSTH; negative thresholds selects the full window
addParamValue(prs,'jitterdefinition','all',...
    @(s)ischar(s)|ismember(s,{'all','burst'}))   % controls the definition of 'jitter'
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
addParamValue(prs,'rasterwindow',[],...
    @(s)isempty(s)|(isnumeric(s)&isequal(length(s),2)))  % time window for raster plot, in seconds
parse(prs,cellid,varargin{:})
g = prs.Results;
dt = 0.0005;   % fixed parameter for now
H = struct('H_psth',{},'H_raster',{});   % initialize structure for figure handles

% Filter events
valid_trials = filterTrials(cellid,'event_type',g.event_type,'event',g.event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);

% Find stimulation peak (latency)
% [lim1 lim2 L at B M] = findstimperiod2(cellid,'event',g.event,...
%     'valid_trials',valid_trials);   %#ok<ASGLU> % find stimulated period and peak time
[~, ~, ~, ~, spt, stats] = ultimate_psth(cellid,g.event_type,g.event,g.window,...
    'dt',dt,'display',g.display,'isadaptive',g.isadaptive,'maxtrialno',Inf,...
    'event_filter',g.event_filter,'filterinput',g.filterinput,...
    'baselinewin',g.baselinewin,'testwin',g.testwin,'relative_threshold',g.relative_threshold,'margin',[-0.01 0.01]);
if g.display
    H(1).H_psth = gcf;   % figure handle for PSTH
end
lim1 = stats.activation_start;  % window start for evoked spikes
lim2 = stats.activation_end;   % window end for evoked spikes
L = stats.activation_peak;   % LATENCY
B = stats.baseline;  % baseline FR
M = stats.maxvalue;  % peak FR

% Evoked spikes
tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],g.event_type,g.event,'valid_trials',valid_trials);   % convert period to epochs relative to event
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find evoked spikes

% Reliability
[R spno selts_evoked1st] = efficiency(cellid,selts_evoked,g.event_type,g.event,'valid_trials',valid_trials);   % efficiency of evoking spikes

% Standard deviation of evoked spikes (jitter)
switch g.jitterdefinition
    case 'all'
        relevokedtimes = abs2reltimes(cellid,selts_evoked,g.event_type,g.event);  % convert absolute spike times to times rel. to the event
    case 'burst'
        relevokedtimes = abs2reltimes(cellid,selts_evoked1st,g.event_type,g.event);  % convert absolute spike times to times rel. to the event
end
J = std(relevokedtimes);   % JITTER

% Plot raster
if g.display
    time = g.window(1):dt:g.window(2);  % time stamps
    tno = size(spt,1);   % number of trials
    H(1).H_raster = figure;
%     pause(.05)
%     jFrame = get(handle(H(1).H_raster),'JavaFrame');   % maximize figure
%     jFrame.setMaximized(true);
    if ~isempty(g.rasterwindow)   % restrict raster window
        rasterinx = time >= g.rasterwindow(1) & time <=g.rasterwindow(2);
        rasterplot(spt(:,rasterinx),time(rasterinx),gcf);  % plot raster
    else
        rasterplot(spt,time,gcf);  % plot raster
    end
    colormap('bone')   % white on black
    set(gcf,'Color','k')
%     set(gcf,'NumberTitle','off','MenuBar','none')
    pause(.05)
    line([lim1 lim1],[1 tno],'LineWidth',1,'Color',[0.8 0 0])
    line([lim2 lim2],[1 tno],'LineWidth',1,'Color',[0.8 0 0])
    pause(.05)
    set(gca,'XColor','w','box','off','TickDir','out')
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'InvertHardcopy','off')
    set(gca,'XColor','w','YColor','w')
    xlabel('Time (s)','Color','w','FontSize',12)
    ylabel('#Trials','Color','w','FontSize',12)
end