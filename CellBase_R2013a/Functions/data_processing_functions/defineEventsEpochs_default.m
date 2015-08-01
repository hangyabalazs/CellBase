function [events,epochs] = defineEventsEpochs_default
%DEFINEEVENTSEPOCHS_DEFAULT   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_DEFAULT defines events and epochs
%   for spike extraction. 
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
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_PULSEON.

%   Edit log: BH 7/6/12

% Define events and epochs
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'PulseOn',       'PulseOn',         'PulseOn',         [-3 3]};    i = i + 1;
events(i,:) = {'PulseOff',      'PulseOff',        'PulseOff',        [-3 3]};    i = i + 1;
events(i,:) = {'BurstOn',       'BurstOn',         'BurstOn',         [-6 6]};    i = i + 1;
events(i,:) = {'BurstOff',      'BurstOff',        'BurstOff',        [-6 6]};    i = i + 1;

% Variable events
events(i,:) = {'PulsePeriod',   'PulseOn',         'PulseOff',        [-1 1]};    i = i + 1;
events(i,:) = {'BurstPeriod',   'BurstOn',         'BurstOff',        [-1 1]};    i = i + 1;
events(i,:) = {'PrePulseIPI',   'PrevPulseOff',    'PulseOn',         [-1 1]};    i = i + 1;
events(i,:) = {'NextPulseIPI',  'PulseOff',        'NextPulseOn',     [-1 1]};    i = i + 1;

% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent     Window             RealWindow
i = 1;

epochs(i,:) = {'FixedBaseline', 'PulseOn',         [-0.005 0.0],       'PrePulseIPI'};    i = i + 1;
epochs(i,:) = {'FixedLightResponse','PulseOn',     [0.0 0.005],        'PulsePeriod'};    i = i + 1;

epochs(i,:) = {'Baseline',       'PrePulseIPI',    [NaN NaN],          'NaN'};            i = i + 1;
epochs(i,:) = {'LightResponse',  'PulsePeriod',    [NaN NaN],          'NaN'};            i = i + 1;