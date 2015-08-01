function [events,epochs] = defineEventsEpochs_laserstim
%DEFINEEVENTSEPOCHS_LASERSTIM   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_LASERSTIM defines events and
%   epochs for spike extraction. 
%
%   EVENTS is a Nx4 cell array with columns corresponding to EventLabel,
%   EventTrigger1, EventTrigger2, Window. EventLabel is the name for
%   referencing the event. EventTrigger1 and EventTrigger2 are names of
%   StimEvent variables (e.g. 'BurstOn'). For fixed windows, the two
%   events are the same; for variable windows, they correspond to the start
%   and end events. Window specifies time offsets relative to the events;
%   e.g. events(1,:) = {'BurstOn', 'BurstOn', 'BurstOn', [-6 6]};
%
%   EPOCH is a Nx4 cell array with columns corresponding to  EpochLabel, 
%   ReferenceEvent, Window, RealWindow. EventLabel is the name for 
%   referencing the epoch. ReferenceEvent should match an EventLabel in 
%   EVENTS (used for calculating the epoch rates). RealWindow is currently
%   not implemented (allocated for later versions).
%
%   DEFINEEVENTSEPOCHS_LASERSTIM defines stimulus events and epochs for
%   photo-stimulation (LASERSTIMPROTOCOL_NI2).
%
%   See also PREALIGNSPIKES and DEFINEEVENTSEPOCHS_DEFAULT.

%   Edit log: BH 7/18/13

% Define events and epochs
%              EventLabel       EventTrigger1    EventTrigger2  Window
i = 1;
events(i,:) = {'BurstOn',       'BurstOn',      'BurstOn',      [-6 6]};   i = i + 1;  
events(i,:) = {'PreBurstIBI',   'PrevBurstOff', 'BurstOn',      [0 0]};    i = i + 1;  
events(i,:) = {'BurstPeriod',   'BurstOn',      'BurstOff',     [0 0]};    i = i + 1;  
events(i,:) = {'NextBurstIBI',  'BurstOff',     'NextBurstOn',  [0 0]};    i = i + 1;  

% Variable events
events(i,:) = {'PreBurstIBI2',  'PrevBurstOff', 'BurstOn',      [-1 0]};   i = i + 1;  
events(i,:) = {'BurstPeriod2',  'BurstOn',      'BurstOff',     [-1 0]};   i = i + 1;  
events(i,:) = {'NextBurstIBI2', 'BurstOff',     'NextBurstOn',  [-1 0]};   i = i + 1;  
events(i,:) = {'OmitPulse',     'OmitPulse',    'OmitPulse',    [-6 6]};   i = i + 1;  
events(i,:) = {'ZeroPulse',     'ZeroPulse',    'ZeroPulse',    [-6 6]};   i = i + 1;

% Define epochs for rate calculations
%               EpochLabel             ReferenceEvent  FixedWindow          RealWindow
i = 1;
epochs(i,:) = {'FixedBaseline5',       'BurstOn',      [-0.005 0.0],  'PreBurstIBI'};   i = i + 1;
epochs(i,:) = {'FixedLightResponse5',  'BurstOn',      [0.0 0.005],   'BurstPeriod'};   i = i + 1;
epochs(i,:) = {'FixedBaseline5a',      'BurstOn',      [-0.005 0.0],  'BurstOn'};       i = i + 1;
epochs(i,:) = {'FixedLightResponse5a', 'BurstOn',      [0.0 0.005],   'BurstOn'};       i = i + 1;
epochs(i,:) = {'FixedBaseline10',      'BurstOn',      [-0.01 0.0],   'PreBurstIBI'};   i = i + 1;
epochs(i,:) = {'FixedLightResponse10', 'BurstOn',      [0.0 0.01],    'BurstPeriod'};   i = i + 1;
epochs(i,:) = {'FixedBaseline20',      'BurstOn',      [-0.02 0.0],   'PreBurstIBI'};   i = i + 1;
epochs(i,:) = {'FixedLightResponse20', 'BurstOn',      [0.0 0.02],    'BurstPeriod'};   i = i + 1;
epochs(i,:) = {'FixedBaseline40',      'BurstOn',      [-0.04 0.0],   'PreBurstIBI'};   i = i + 1;
epochs(i,:) = {'FixedLightResponse40', 'BurstOn',      [0.0 0.04],    'BurstPeriod'};   i = i + 1;
epochs(i,:) = {'BurstBaseline',        'PreBurstIBI',  [NaN NaN],     'NaN'};           i = i + 1;
epochs(i,:) = {'BurstResponse',        'BurstPeriod',  [NaN NaN],     'NaN'};           i = i + 1;