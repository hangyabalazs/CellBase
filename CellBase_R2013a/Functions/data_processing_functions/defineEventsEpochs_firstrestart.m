function [events,epochs] = defineEventsEpochs_firstrestart
%DEFINEEVENTSEPOCHS_FIRSTRESTART   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_FIRSTRESTART defines events and
%   epochs for spike extraction. 
%
%   EVENTS is a Nx4 cell array with columns corresponding to EventLabel,
%   EventTrigger1, EventTrigger2, Window. EventLabel is the name for
%   referencing the event. EventTrigger1 and EventTrigger2 are names of
%   TrialEvent variables (e.g. 'LeftPortIn'). For fixed windows, the two
%   events are the same; for variable windows, they correspond to the start
%   and end events. Window specifies time offsets relative to the events;
%   e.g. events(1,:) = {'OdorValveOn','OdorValveOn','OdorValveOn',[-3 3]};
%
%   EPOCH is a Nx4 cell array with columns corresponding to  EpochLabel, 
%   ReferenceEvent, Window, RealWindow. EventLabel is the name for 
%   referencing the epoch. ReferenceEvent should match an EventLabel in 
%   EVENTS (used for calculating the epoch rates). RealWindow is currently
%   not implemented (allocated for later versions).
%
%   DEFINEEVENTSEPOCHS_FIRSTRESTART defines an event for the first
%   ITI-restarting lick.
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_DEFAULT.

%   Edit log: BH 9/18/13

% Define events and epochs - variable events
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'FirstRestart',  'FirstRestart',   'FirstRestart',     [-6 6]};
epochs = [];