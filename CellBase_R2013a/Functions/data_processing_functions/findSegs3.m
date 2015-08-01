function tseg = findSegs3(cellid,varargin)
%FINDSEGS3   Find segments of a session.
%   TSEG = FINDSEGS3(CELLID) finds time segments for a given cell (CELLID)
%   and returns the start (first row of TSEG) and end timestamps (second
%   row of TSEG); thus, every column of TSEG corresponds to the start and
%   end of a segment. Segments are selected based on a particular filter,
%   which can be set via a 'segfilter', filter type parameter value pair:
%   TSEG = FINDSEGS3(CELLID,'segfilter',FILTERNAME). The default FILTERNAME
%   is 'stim_excl' (see below).
%
%   Filters currently implemented (possible values for FILTERNAME):
%   	'stim_excl' - gaps between protocols excluding pulses during
%           behavior
%     	'prestim' - pre-stimulus segments: from the offset of a pulse to
%           the onset of the next
%     	'prestim2' - pre-stimulus segments: fixed window before each
%           stimulus onset; window is set by 'prepulseinterval' input
%           agrument (default, 1s); note that 'margin', 'min_int' and 
%           'max_int' arguments are disregarded
%     	'stim' - stimulation segments during stimulation and behavior
%           protocols; note that 'margin', 'min_int' and 'max_int'
%           arguments are disregarded
%     	'spont' - gaps between behav/stim protocols
%       'stim_excl_no_behavior' - segments between protocols and between
%           light pulses
%       'stim_excl_nb' - windows around light pulses both during and
%           outside behavior epochs are excluded; the windows are defined 
%           by the 'light_stimulation_duration' variable
%       'stimfb_excl_nb' - windows around light pulses both during and
%           outside behavior epochs are excluded; the windows are defined 
%           by the 'light_stimulation_duration' variable; also, windows
%           around behavioral feedback (reward or punishment) are excluded;
%           the windows are defined by the 'feedback_duration' variable
%       'fb_incl_nb' - only windows around behavioral feedback (reward or
%           punishment) are included; the windows are defined by the
%           'feedback_duration' variable
%   To implement a new filter, a subfunction in the form [TSEG G] =
%   FILTERNAME(G,CELLID) should be appended to the code. G is the structure
%   containing the input arguments to FINDSEGS3. Returning G from the
%   subfunction makes it possible to force fixed values onto these
%   arguments in a filter-dependent manner.
%
%   Additional optional input arguments to FINDSEGS3 (parameter, value
%   pairs, with their default values):
%   	'margins', [1 -1] - define a margin around the events that will be
%           excluded; in seconds
%       'min_int', 1 - minimum interval length of time segment; in seconds
%       'max_int', Inf - maximum interval length of time segment; in seconds
%       'light_activation_duration', [0 1] - time window after light 
%           activation; in seconds
%       'prepulseinterval', 1 - % time window to include before
%           stimulation; in seconds (for prestim2 filter)
%
%   See also EXTRACTSEGSPIKES.

%   Edit log: SPR 4/10, BH 4/25/12, 7/5/13

% Default arguments
default_args = {...
    'margins',      [1 -1];...       % in seconds
    'min_int',      1;...            % minimum interval of time segment (secs)
    'max_int',      Inf;...          % max interval
    'segfilter',    'stim_excl';...  % default: exclude stimulation
    'light_activation_duration',   [0 1];...   % time window after stimulation, in seconds
    'feedback_duration',   [-0.1 0.1];...   % time window around feedback, in seconds (for stimfb_excl_nb filter)
    'prepulseinterval',   1;...      % time window to include before stimulation, in seconds (for prestim2 filter)
    };
[g,error] = parse_args(default_args,varargin{:});

% Apply filter
[tseg g] = feval(str2func(g.segfilter),g,cellid);

% Filter according to margin and interval
tseg = tseg + repmat(g.margins,length(tseg),1)';
tseg = tseg(:,~isnan(diff(tseg))&diff(tseg)>g.min_int&diff(tseg)<g.max_int);

% -------------------------------------------------------------------------
function [all_segs g] = stim_excl_no_behavior(g,cellid) %#ok<*DEFNU>

% Major gaps between protocols
major_segs = spont(g,cellid);

% Pre-stimulus segments: from the offset of a pulse to the onset of the next
prestim_segs = prestim(g,cellid);

% Segments between protocols and between light pulses
all_segs = [major_segs prestim_segs];

% -------------------------------------------------------------------------
function [prestim_segs g] = prestim(g,cellid)

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');

    % Pre-stimulus segments: from the offset of a pulse to the onset of the next
    if ~isempty(SE.PulseOn)
        prestim_segs = [SE.PrevPulseOff; SE.PulseOn];
    else
        prestim_segs = [];
    end
catch
    disp([cellid ': There was no stim protocol for this session.'])
    prestim_segs = [];
end

% -------------------------------------------------------------------------
function [prestim_segs g] = prestim2(g,cellid)

% Change arguments
g.margins = [0 0];
g.min_int = 0;
g.max_int = Inf;

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');

    % Pre-stimulus segments: from the offset of a pulse to the onset of the next
    if ~isempty(SE.PulseOn)
        prestim_segs = [SE.PulseOn-g.prepulseinterval; SE.PulseOn];
    else
        prestim_segs = [];
    end
catch
    disp([cellid ': There was no stim protocol for this session.'])
    prestim_segs = [];
end

% -------------------------------------------------------------------------
function [prestim_segs g] = prestim3(g,cellid)

% Change arguments
g.margins = [0 0];
g.min_int = 0;
g.max_int = Inf;

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');

    % Pre-stimulus segments: from the offset of a pulse to the onset of the next
    if ~isempty(SE.BurstOn)
        prestim_segs = [SE.BurstOn-g.prepulseinterval; SE.BurstOn];
    else
        prestim_segs = [];
    end
catch
    disp([cellid ': There was no stim protocol for this session.'])
    prestim_segs = [];
end

% -------------------------------------------------------------------------
function [major_segs g] = spont(g,cellid)

% Load events
Ev = loadcb(cellid,'Events');

% Session begin and end
sess_begin = Ev.Events_TimeStamps(1);
sess_end = Ev.Events_TimeStamps(end);

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');
    
    % Start and end points of stimulation protocols 
    stim_begin = unique(SE.ProtocolStart(~isnan(SE.ProtocolStart)));
    stim_end = unique(SE.ProtocolEnd(~isnan(SE.ProtocolEnd)));
    
    % No protocol end event
    if isempty(stim_end)
        stim_end = SE.PulseOn(end) + 2;   % 2 seconds after the last pulse
    end
catch
    disp([cellid ': There was no stim protocol for this session.'])
    stim_begin = [];
    stim_end = [];
end

% Trial events
try
    % Load trial events
    TE = loadcb(cellid,'TrialEvents');
    
    % Start and end points of behavior protocols
    behav_begin = unique(TE.ProtocolStart(~isnan(TE.ProtocolStart)));
    behav_end = unique(TE.ProtocolEnd(~isnan(TE.ProtocolEnd)));
    if length(behav_begin) > 1
        error('FindSegs does not work for sessions with multiple behavior protocols.')
    end
catch
    disp([cellid ': There was no behavioral protocol for this session.'])
    behav_begin = [];
    behav_end = [];
end

% Start and end points of stimulus and behavior protocols
prot_begin = [stim_begin behav_begin];
prot_end = [stim_end behav_end];

% Major gaps between protocols
major_segs = sort([sess_begin sess_end prot_begin prot_end]);
major_segs = reshape(major_segs,2,[]);

% -------------------------------------------------------------------------
function [tseg g] = stim(g,cellid)

% Change arguments
g.margins = [0 0];
g.min_int = 0;
g.max_int = Inf;

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');
    
    % Druing stimulation protocols: from pulse onset, using light activation duration argument
    stim_segs = [SE.PulseOn + g.light_activation_duration(1); ...
        SE.PulseOn + g.light_activation_duration(2)];
catch
    disp([cellid ': There was no stim protocol for this session.'])
    stim_segs = [];
end

% Behav. segments, excluding pulses during behavior
[all_behav_segs all_behav_stimsegs] = allbehavsegs(g,cellid);

% Stimulation segments during stimulation and behavior protocols
tseg = [stim_segs all_behav_stimsegs];

% -------------------------------------------------------------------------
function [all_behav_segs all_behav_stimsegs] = allbehavsegs(g,cellid)

% Load trial events
try
    % Load trial events
    TE = loadcb(cellid,'TrialEvents');
    
    % Exclude pulses during behavior
    all_behav_stims = unique([TE.AllPulseOns{~isempty(TE.AllPulseOns)} TE.PulseOn(~isnan(TE.PulseOn))]);
    all_behav_stimsegs = [all_behav_stims + g.light_activation_duration(1); ...
        all_behav_stims + g.light_activation_duration(2)];  % stimulus segments
    all_behav_segs = reshape(sort([TE.ProtocolStart(1) all_behav_stimsegs(:)' TE.ProtocolEnd(1)]),2,[]);  % behav. segments excluding stimulus segments
catch
    disp([cellid ': There was no behavioral protocol for this session.'])
    all_behav_segs = [];
    all_behav_stimsegs = [];
end

% -------------------------------------------------------------------------
function [tseg g] = stim_excl(g,cellid)

% Major gaps between protocols
major_segs = spont(g,cellid);

% Behav. segments, excluding pulses during behavior
all_behav_segs = allbehavsegs(g,cellid);

% Gaps between protocols excluding pulses during behavior
tseg = sort([major_segs all_behav_segs],2);

% -------------------------------------------------------------------------
function [tseg g] = stim_excl_nb(g,cellid)

% Load events
Ev = loadcb(cellid,'Events');

% Session begin and end
sess_begin = Ev.Events_TimeStamps(1);
sess_end = Ev.Events_TimeStamps(end);

% Stimulus events
try
    % Load stimulus events
    SE = loadcb(cellid,'StimEvents');
    
    % Druing stimulation protocols: from pulse onset, using light activation duration argument
    stim_segs = [SE.PulseOn + g.light_activation_duration(1); ...
        SE.PulseOn + g.light_activation_duration(2)];
catch
    disp([cellid ': There was no stim protocol for this session.'])
    stim_segs = [];
end

% Light stimulation during behavior
try
    % Load trial events
    TE = loadcb(cellid,'TrialEvents');
    
    % Druing stimulation protocols: from pulse onset, using light activation duration argument
    TE.LightStimulation2 = TE.LightStimulation2(~isnan(zero2nan(TE.LightStimulation2)));
    stim_segs_behav = [TE.LightStimulation2 + g.light_activation_duration(1); ...
        TE.LightStimulation2 + g.light_activation_duration(2)];
    if ~isempty(stim_segs_behav)
        keyboard   % this needs debugging
    end
    
catch
    disp([cellid ': There was no behavioral protocol for this session.'])
    stim_segs_behav = [];
end

% Stimulation segments during stimulation and behavior protocols
all_stim_segs = [stim_segs stim_segs_behav];

% Gaps between stim. segments
tseg = [sess_begin all_stim_segs(2,:); all_stim_segs(1,:) sess_end];

% -------------------------------------------------------------------------
function [tseg g] = stimfb_excl_nb(g,cellid)

% Exclude stimulation
[tseg g] = stim_excl_nb(g,cellid);   % segments between stimulation windows
sess_begin = tseg(1,1);  % extract the session limits: session begin and end
sess_end = tseg(2,end);
stim_segs = [tseg(2,1:end-1); tseg(1,2:end)];   % extract the stimulation windows

% Checking whether 'DeliverFeedback' event is available
sesstype = getvalue('session_type',cellid);
if isequal(sesstype,{'feedbackdelay'})
    alignevent_fa = 'DeliverFeedback';
else
    alignevent_fa = 'LeftPortIn';
end
if isequal(sesstype,{'feedbackdelay'})
    alignevent_hit = 'DeliverFeedback';
else
    alignevent_hit = 'LeftWaterValveOn';
end

% Load trial events
try
    % Load trial events
    TE = loadcb(cellid,'TrialEvents');
    
    % Segments to exclude
    hit_trials = ~isnan(TE.Hit);  % hit trials
    exclude_time_hit = TE.TrialStart(hit_trials) + TE.(alignevent_hit)(hit_trials);
    fa_trials = ~isnan(TE.FalseAlarm);  % false alarm trials
    exclude_time_fa = TE.TrialStart(fa_trials) + TE.(alignevent_fa)(fa_trials);
    exclude_time = sort([exclude_time_hit exclude_time_fa]);
    behav_segs = [exclude_time + g.feedback_duration(1); ...
        exclude_time + g.feedback_duration(2)];
catch
    disp([cellid ': There was no behavioral protocol for this session.'])
    behav_segs = [];
end

% Final segments
allseg = [stim_segs behav_segs];
allseg = sortrows(allseg',1)';
tseg = [sess_begin allseg(2,:); allseg(1,:) sess_end];

% -------------------------------------------------------------------------
function [behav_segs g] = fb_incl_nb(g,cellid)

% Checking whether 'DeliverFeedback' event is available
sesstype = getvalue('session_type',cellid);
if isequal(sesstype,{'feedbackdelay'})
    alignevent_fa = 'DeliverFeedback';
else
    alignevent_fa = 'LeftPortIn';
end
if isequal(sesstype,{'feedbackdelay'})
    alignevent_hit = 'DeliverFeedback';
else
    alignevent_hit = 'LeftWaterValveOn';
end

% Load trial events
try
    % Load trial events
    TE = loadcb(cellid,'TrialEvents');
    
    % Segments to include
    hit_trials = ~isnan(TE.Hit);  % hit trials
    time_hit = TE.TrialStart(hit_trials) + TE.(alignevent_hit)(hit_trials);
    fa_trials = ~isnan(TE.FalseAlarm);  % false alarm trials
    time_fa = TE.TrialStart(fa_trials) + TE.(alignevent_fa)(fa_trials);
    time_fb = sort([time_hit time_fa]);
    behav_segs = [time_fb + g.feedback_duration(1); ...
        time_fb + g.feedback_duration(2)];
catch
    disp([cellid ': There was no behavioral protocol for this session.'])
    behav_segs = [];
end

% -------------------------------------------------------------------------
function [behav_segs g] = cue_incl_nb(g,cellid)

% align to cue onset
alignevent = 'StimulusOn';

% Load trial events
try
    % Load trial events
    TE = loadcb(cellid,'TrialEvents');
    
    % Segments to include
    time_pc = TE.TrialStart + TE.(alignevent);
    behav_segs = [time_pc + g.feedback_duration(1); ...
        time_pc + g.feedback_duration(2)];
catch
    disp([cellid ': There was no behavioral protocol for this session.'])
    behav_segs = [];
end