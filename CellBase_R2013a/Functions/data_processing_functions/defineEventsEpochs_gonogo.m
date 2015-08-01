function [events,epochs] = defineEventsEpochs_gonogo
%DEFINEEVENTSEPOCHS_GONOGO   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_GONOGO defines events and epochs
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
%   DEFINEEVENTSEPOCHS_GONOGO defines events and epochs for auditory
%   go-nogo task.
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_DEFAULT.

%   Edit log: BH 7/6/12

% Define events and epochs
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'StimulusOn',    'StimulusOn',      'StimulusOn',      [-6 6]};    i = i + 1;
events(i,:) = {'StimulusOff',   'StimulusOff',     'StimulusOff',     [-6 6]};    i = i + 1;
events(i,:) = {'LeftPortIn',    'LeftPortIn',      'LeftPortIn',      [-6 6]};    i = i + 1;
events(i,:) = {'LeftWaterValveOn','LeftWaterValveOn','LeftWaterValveOn',[-6 6]};  i = i + 1;
events(i,:) = {'DeliverFeedback','DeliverFeedback','DeliverFeedback', [-6 6]};    i = i + 1;
events(i,:) = {'LickIn',        'LickIn',          'LickIn',          [-6 6]};    i = i + 1;
events(i,:) = {'LickOut',       'LickOut',         'LickOut',         [-6 6]};    i = i + 1;

% Variable events
events(i,:) = {'StimulusSampling','StimulusOn',    'StimulusOff',     [-6 6]};    i = i + 1;

% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent      FixedWindow       RealWindow
i = 1;
epochs(i,:) = {'StimulusOn',    'StimulusOn',       [0.0 0.3],        'StimulusSampling'};    i = i + 1;