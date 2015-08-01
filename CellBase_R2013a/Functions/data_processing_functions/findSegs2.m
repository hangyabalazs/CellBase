function tseg = findSegs2(cellid,varargin)
%FINDSEGS2   Find segments of a session.
%   TSEG = FINDSEGS2(CELLID) finds time segments for a given cell (CELLID)
%   and returns the start (first row of TSEG) and end timestamps (second
%   row of TSEG); thus, every column of TSEG corresponds to the start and
%   end of a segment. Segments are selected based on a particular filter,
%   which can be set via a 'segfilter', filter type parameter value pair:
%   TSEG = FINDSEGS2(CELLID,'segfilter',FILTERNAME). The default FILTERNAME
%   is 'stim_excl' (see below).
%
%   Filters currently implemented (possible values for FILTERNAME):
%   	'stim_excl' - gaps between protocols excluding pulses during
%           behavior
%     	'prestim' - pre-stimulus segments: from the offset of a pulse to
%           the onset of the next
%     	'stim' - stimulation segments during stimulation and behavior
%           protocols; note that 'margin', 'min_int' and 'max_int'
%           arguments are disregarded
%     	'spont' - gaps between behav/stim protocols
%       'stim_excl_no_behavior' - segments between protocols and between
%           light pulses
%   To implement a new filter, a subfunction in the form [TSEG G] =
%   FILTERNAME(G,CELLID) should be appended to the code. G is the structure
%   containing the input arguments to FINDSEGS2. Returning G from the
%   subfunction makes it possible to force fixed values onto these
%   arguments in a filter-dependent manner.
%
%   Additional optional input arguments to FINDSEGS2 (parameter, value
%   pairs, with their default values):
%   	'margins', [1 -1] - define a margin around the events that will be
%           excluded; in seconds
%       'min_int', 1 - minimum interval length of time segment; in seconds
%       'max_int', Inf - maximum interval length of time segment; in seconds
%       'light_activation_duration', [0 1] - time window after light 
%           activation; in seconds
%
%   See also EXTRACTSEGSPIKES.

%   Edit log: SPR 4/10, BH 4/25/12

% Default arguments
default_args = {...
    'margins',      [1 -1];...       % in seconds
    'min_int',      1;...            % minimum interval of time segment (secs)
    'max_int',      Inf;...          % max interval
    'segfilter',    'stim_excl';...  % default: exclude stimulation
    'light_activation_duration',   [0 1];...   % time window after stimulation, in seconds
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

% Load stimulus events
SE = loadcb(cellid,'StimEvents');

% Pre-stimulus segments: from the offset of a pulse to the onset of the next
if ~isempty(SE.PulseOn)
    prestim_segs = [SE.PrevPulseOff; SE.PulseOn];
else
    prestim_segs = [];
end

% -------------------------------------------------------------------------
function [major_segs g] = spont(g,cellid)

% Load events
Ev = loadcb(cellid,'Events');

% Load stimulus events
try
    SE = loadcb(cellid,'StimEvents');
catch
    disp('There was no stim protocol for this session.')
end

% Load trial events
try
    TE = loadcb(cellid,'TrialEvents');
catch
    disp('There was no behavioral protocol for this session.')
end

% Session begin and end
sess_begin = Ev.Events_TimeStamps(1);
sess_end = Ev.Events_TimeStamps(end);

% Start and end points of behavior protocols
behav_begin = unique(TE.ProtocolStart(~isnan(TE.ProtocolStart)));
behav_end = unique(TE.ProtocolEnd(~isnan(TE.ProtocolEnd)));
if length(behav_begin) > 1
    error('FindSegs does not work for sessions with multiple behavior protocols.')
end

% Start and end points of stimulation protocols 
stim_begin = unique(SE.ProtocolStart(~isnan(SE.ProtocolStart)));
stim_end = unique(SE.ProtocolEnd(~isnan(SE.ProtocolEnd)));

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

% Load stimulus events
try
    SE = loadcb(cellid,'StimEvents');
catch
    disp('There was no stim protocol for this session.')
end

% Druing stimulation protocols: from pulse onset, using light activation duration argument
stim_segs = [SE.PulseOn + g.light_activation_duration(1); ...
    SE.PulseOn + g.light_activation_duration(2)];

% Pulses during behavior
[all_behav_segs all_behav_stimsegs] = allbehavsegs(g,cellid);

% Stimulation segments during stimulation and behavior protocols
tseg = [stim_segs all_behav_stimsegs];

% -------------------------------------------------------------------------
function [all_behav_segs all_behav_stimsegs] = allbehavsegs(g,cellid)

% Load trial events
try
    TE = loadcb(cellid,'TrialEvents');
catch
    disp('There was no behavioral protocol for this session.')
end

% Pulses during behavior
all_behav_stims = unique([TE.AllPulseOns{~isempty(TE.AllPulseOns)} TE.PulseOn(~isnan(TE.PulseOn))]); 
all_behav_stimsegs = [all_behav_stims + g.light_activation_duration(1); ...
    all_behav_stims + g.light_activation_duration(2)];  % stimulus segments
all_behav_segs = reshape(sort([TE.ProtocolStart(1) all_behav_stimsegs(:)' TE.ProtocolEnd(1)]),2,[]);  % behav. segments excluding stimulus segments

% -------------------------------------------------------------------------
function [tseg g] = stim_excl(g,cellid)

% Major gaps between protocols
major_segs = spont(g,cellid);

% Pulses during behavior
all_behav_segs = allbehavsegs(g,cellid);

% Gaps between protocols excluding pulses during behavior
tseg = sort([major_segs all_behav_segs],2);