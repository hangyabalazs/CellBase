function abstimes = rel2abstimes(cellid,reltimes,event_type,event,varargin)
%REL2ABSTIMES   Calculate absolute time windows relative to events.
%   ABSTIMES = REL2ABSTIMES(CELLID,RELTIMES,EVENT_TYPE,EVENT) calculates
%   time stamps (ABSTIMES) corresponding to RELTIMES (time values relative
%   to an EVENT, in seconds) for a given cell (CELLID). EVENT_TYPE should
%   determine whether EVENT is a field of stimulus events ('stim') or trial
%   events ('trial'). Alternatively, RELTIMES can contain time stamps
%   relative to all trials separately in a cell arrray (STIMSPIKES
%   structure, see PREALIGNSPIKES).
%
%   ABSTIMES = REL2ABSTIMES(CELLID,RELTIMES,EVENT_TYPE,EVENT,'VALID_TRIALS',VT)
%   filters the trials according to the 'VALID_TRIALS', TR parameter value
%   pair.
%
%   See also ABS2RELTIMES and FINDSEGS3.

%   Event log: BH 4/25/12, 5/4/12

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'reltimes',@(s)isnumeric(s)|iscell(s))  % time stamps to convert
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@ischar)   % reference event
addParamValue(prs,'valid_trials','all',@(s)isnumeric(s)||islogical(s))  % valid trials - use all trials by default
parse(prs,cellid,reltimes,event_type,event,varargin{:})
g = prs.Results;

% Load event structure
event_type = lower(event_type(1:4));
switch event_type
    case 'stim'
        
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');
        catch ME
            disp('There was no stim protocol for ths session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');
        catch ME
            disp('There was no behavioral protocol for ths session.')
            error(ME.message)
        end
        
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

% Valid trials
valid_trials = parseValidTrials(VE,event,g.valid_trials);

% Calculate time stamps
if isnumeric(reltimes)  % align a set of relative times to all trials
    tno = length(reltimes);
    abstimes = nan(tno,length(find(valid_trials)));     % works for both logical and numeric 'valid_trials'
    if isequal(event_type,'tria') && ~isequal(event,'TrialStart')
        for k = 1:tno
            abstimes(k,:) = VE.(event)(valid_trials) + VE.TrialStart(valid_trials) + reltimes(k);     % trial events other than TrialStart are stored relative to TrialStart
        end
    else
        for k = 1:tno
            abstimes(k,:) = VE.(event)(valid_trials) + reltimes(k);
        end
    end
elseif iscell(reltimes)  % convert prealigned time stamps to absolute times
    reltimes = reltimes(valid_trials);  % restrict to valid trials
    VEev = VE.(event)(valid_trials);
    VEts = VE.TrialStart(valid_trials);
    tno = length(reltimes);
    abstimes = cell(tno,1);
    if isequal(event_type,'tria') && ~isequal(event,'TrialStart')
        for k = 1:tno
            abstimes{k} = VEev(k) + VEts(k) + reltimes{k};     % trial events other than TrialStart are stored relative to TrialStart
        end
    else
        for k = 1:tno
            abstimes{k} = VEev(k) + reltimes{k};
        end
    end
    abstimes = cell2mat(abstimes);
end