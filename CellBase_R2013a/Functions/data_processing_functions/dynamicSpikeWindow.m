function [stimes2 starttimes endtimes] = dynamicSpikeWindow(stimes,VE,TriggerEvent,FirstEvent,LastEvent)
%DYNAMICSPIKEWINDOW   Restrict spike times to event windows.
%   ST = DYNAMICSPIKEWINDOW(ST,VE,TRIGGEREVENT,FIRSTEVENT,LASTEVENT)
%   excludes pre-aligned spikes (ST) before a preceding (FIRSTEVENT) and
%   after a following (LASTEVENT) event. Trial or stimulus events (VE) and
%   the name of the trigger event for the pre-alignment (TRIGGEREVENT)
%   should be passed onto the function. No restriction is applied for empty
%   FIRSTEVENT or LASTEVENT.
%
%   [ST T1 T2] = DYNAMICSPIKEWINDOW(ST,VE,TRIGGEREVENT,FIRSTEVENT,LASTEVENT)
%   also returns start and end time stamps for all trials.
%
%   See also PREALIGNSPIKES.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   12-Dec-2012

%   Edit log: BH 12/12/13

% Window limits
NumTrials = length(stimes);
if ~isempty(LastEvent)   % end times
    endtimes = VE.(LastEvent) - VE.(TriggerEvent);
else
    endtimes = ones(1,NumTrials) * Inf;
end
if ~isempty(FirstEvent)   % start times
    starttimes = VE.(FirstEvent) - VE.(TriggerEvent);
else
    starttimes = ones(1,NumTrials) * -Inf;
end

% Restrict spike times
stimes2 = cell(size(stimes));
for iT = 1:NumTrials
    lst = stimes{iT};
    stimes2{iT} = lst(lst>starttimes(iT)&lst<endtimes(iT));
end