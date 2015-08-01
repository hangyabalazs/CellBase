function EventTimes = trialevents2relativetime(TE,RefEvent,Events)
%TRIALEVENTS2RELATIVETIME   Calculate time values relative to reference event.
%   EVENTTIMES = TRIALEVENTS2RELATIVETIME(TE,REFEVENT,EVENTS) calculates
%   relative times (EVENTTIMES) of trigger events (EVENTS) with respect to
%   reference events (REFEVENT). Input arguments are field names of the
%   trial events structure (TE). EVENTS can be a cell array of strings as
%   well.
%
%   See also VIEWCELL2B.

%   Edit log: BH 6/23/11

% Input arguments
if ischar(Events)
    Events = {Events};  % only one event
end

% Calculate relative times
NUMtrials = length(TE.(Events{1}));
NUMevents = length(Events);
if ~iscell(TE.(RefEvent))
    if ~any(cellfun(@iscell,arrayfun(@(s)TE.(Events{s}),1:NUMevents,'UniformOutput',false)))
        EventTimes = nan(NUMevents,NUMtrials);
        for iE = 1:length(Events)
            TriggerEvent = Events{iE};
            EventTimes(iE,:) = TE.(TriggerEvent) - TE.(RefEvent);
        end
    else   % e.g. if showEvents is LickIn
        EventTimes = cell(NUMevents,NUMtrials);
        for iE = 1:length(Events)
            TriggerEvent = Events{iE};
            for iT = 1:NUMtrials
                if iscell(TE.(TriggerEvent))
                    EventTimes{iE,iT} = TE.(TriggerEvent){iT} - TE.(RefEvent)(iT);
                else
                    EventTimes{iE,iT} = TE.(TriggerEvent)(iT) - TE.(RefEvent)(iT);
                end
            end
        end
    end    
else   % deals with multiple events per trial, e.g. lick-aligned raster - BH
    for iE = 1:length(Events)
        TriggerEvent = Events{iE}; 
        EventTimes = cell(1,NUMevents);
        if ~iscell(TE.(TriggerEvent))
            EventTimes{iE} = arrayfun(@(s)TE.(TriggerEvent)(s)-TE.(RefEvent){s},(1:NUMtrials)','UniformOutput',false);
        else    % e.g. if showEvents is LickOut (problem: sometimes LickOut is in the following trial -> numbers disagree)
%             EventTimes{iE} = arrayfun(@(s)TE.(TriggerEvent){s}-TE.(RefEvent){s},(1:NUMtrials)','UniformOutput',false);
            error('ShowEvents for mutiple TriggerEvents and multiple ShowEvents at the same time has not been implemented.')
        end
    end
end