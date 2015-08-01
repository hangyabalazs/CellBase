function run_makeapsth(cellid,varargin)
% FUNCTION RUN_MAKEAPSTH(CELLID)
% Makes a PSTH with adaptive Gaussian kernel for HomeZoneOut1.

% For an event,
% Get the aligned spikes for the event.
% Get the Previous and Next event windows.
% Make binraster.
% Make trial partitions (i.e. conditions for which we want separate PSTHs).
% MAke APSTH for each partition.
% Plot.


if nargin<1,
    cellid = 'd046_100816a_3.2';
end

default_args={...
    'window',               [-2 3];...
    'dt',                   0.001;...
    'sigma',                0.04;...
    'FigureNum',            1;...
    'TriggerName',          'NextForwardRun3';...
    'SortEvent',            'NextTriggerZoneIn';...
    'LastEvents',           '';...
    'eventtype',             'behav';... % 'behav'
    'ShowEvents',           {{'NextTriggerZoneIn'}};...
    'ShowEventsColors',     {{'r'}};...
    'Num2Plot',             'all';...
    ...
    'PlotDashedEvent',      '';...
    'PlotDashedCondition',  'min';...
    'PSTHPlot',             1;...
    'PSTHlinewidth',        1.5;...
    'DashedLineStyle',  ':';...
    'Partitions',           '#SoundID';...
    'PrintCellID',          'on';...
    'PrintCellIDPos',       'bottom-right';...
    };
[g,error]=parse_args(default_args,varargin{:});

margin = g.sigma*3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array

% load
SP = loadcb(cellid,'Eventspikes');
TE = loadcb(cellid,'TrialEvent');

% Find position of TriggerName.
trigger_pos = findcellstr(SP.events(:,1),g.TriggerName);

% TriggerName mismatch
if (trigger_pos == 0)
    error('Trigger name not found');
else
    % Get TriggerEvent which corresponds to a field in TrialEvents.
    TriggerEvent=SP.events{trigger_pos,2};
end

if ~isfield(TE,TriggerEvent),
    error('TriggerEvent mismatch: supply correct Events structure')
end

alltrials = 1:size(SP.event_stimes{1},2);
stimes  = SP.event_stimes{trigger_pos}(alltrials);
% windows = SP.event_windows{trigger_pos}(:,alltrials);

if ~iscellstr(g.LastEvents) && (strcmpi(g.LastEvents,'none') || isempty(g.LastEvents))
    window_margin = SP.events{trigger_pos,4};
    ev_windows = SP.event_windows{trigger_pos};
else
    window_margin = [g.window(1)-2*g.dt 0];
    ev_windows = get_last_evtime(TE,TriggerEvent,g.LastEvents);
end

%%% MAKE THE MAIN RASTER
binraster = stimes2binraster(stimes,time,g.dt,ev_windows,window_margin);
figure
imagesc(binraster);

[COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
subplot(211)
imagesc(binraster(COMPTRIALS{1},:))
subplot(212)
imagesc(binraster(COMPTRIALS{2},:))
valid_trials=find(~isnan(getfield(TE,TriggerEvent)));

% valid_trials=valid_trials(0.1<TE.NextPulseOn(valid_trials)-TE.PulseOn(valid_trials)>0.2);
%%% Could be put as an option
% if g.OdorPairID == 0
%     valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.ValidTrials));
% else
%     valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.OdorPairID,g.ValidTrials));
% end

% [psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);

[psth, spsth, spsth_se] = binraster2apsth(binraster,g.dt,time,g.sigma,COMPTRIALS,valid_trials);

figure(g.FigureNum)
plot(time,psth)