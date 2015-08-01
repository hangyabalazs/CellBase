function ev_window = get_last_evtime(TE,TriggerEvent,LastEvent)
%GET_LAST_EVTIME   Event windows.
%   EV_WINDOW = GET_LAST_EVTIME(TE,TRIGGEREVENT,LASTEVENT) creates event
%   windows using times between TRIGGEREVENTS and LASTEVENTS. Trial event
%   structure should be provided in TE.
%
%   See also VIEWCELL2B.

%   Edit log: AK 07/1, BH 6/23/11

% Convert to cell array of strings
if iscellstr(LastEvent)
    N = length(LastEvent);
else
    LastEvent = cellstr(LastEvent);
    N = 1;
end

% LastEvent times relative to TriggerEvent
for iE = 1:N
    Ltime(iE,:) = TE.(LastEvent{iE}) - TE.(TriggerEvent);
end

% Create 'ev_window'
ev_window = zeros(length(Ltime),2);
if N == 1
    ev_window(:,2) = Ltime';    
else
    ev_window(:,2) = min(Ltime)';
end