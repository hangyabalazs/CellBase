function  prealignSpikes(mycellids,varargin)
%PREALIGNSPIKES   Align spikes to events and calculates rates for each epoch.
%   PREALIGNSPIKES(MYCELLIDS) alignes spikes of cells corresponding to 
%   MYCELLIDS to multiple events. To results are saved in a .mat file (STIMES).
%   STIMES data format:
%   ST.CELLID - cell ID
%   ST.EVENTS - EventLabel, EventTrigger1, EventTrigger2, Margin (see
%       DEFINEEVENTSEPOCHS_DEFAULT)
%   ST.EVENT_STIMES - spike times for corresponding events
%   ST.EVENT_WINDOWS - time windows for corresponding events
%   ST.EPOCHS - EpochLabel, ReferenceEvent, FixedWindow, RealWindow (see
%       DEFINEEVENTSEPOCHS_DEFAULT)
%   ST.EPOCH_RATES - firing rates for corresponding epochs
%   If MYCELLIDS is 'all', then the function will run on all the cells in
%   CellBase.
%
%   Additional parameter value pairs can be passed to the function.
%   Optional input arguments (with default values):
%       'ifsave', true - control saving
%       'ifappend', true - if true, append new events to preexisting file
%       'writing_behavior', [] - if not empty, take precedence over 'ifappend' 
%           'append' - append new events to preexisting file
%           'replace' - attempt to replace new events only
%           'overwrite' - replace file
%       'FUNdefineEventsEpochs', @defineEventsEpochs_default - function that
%           is called to define events (you may create your own
%           event/epoch def. file, eg.: @defineEventsEpochs_adam1106)
%       'filetype', 'event' - ('behav') for behavior and 'stim' for light
%           stimulation
%       'events', [] - events can be passed directly; definition file is
%           not used
%       'epochs', [] - epochs can be passed directly; definition file is
%           not used
%
%   Examples:
%   prealignSpikes(cellid,'FUNdefineEventsEpochs',...
%       @defineEventsEpochs_gonogo,'filetype','event','ifsave',1,'ifappend',0)
%   prealignSpikes(CELLIDLIST(190),'events',...
%       {'PulseOnS' 'PulseOn' 'PulseOn' [0.0086 0.01]},'epochs',[],...
%       'filetype','stim','ifsave',1,'ifappend',1)
%
%   See also VIEWCELL2B and PLOT_RASTER2A.

%   Edit log: AK 11/06; 2/11, BH 6/27/11, 5/10/12

% Input arguments
prs = inputParser;
addRequired(prs,'mycellids',@(s)iscell(s)|ischar(s))    % cell array of cell IDs
addParamValue(prs,'ifsave',false,@(s)islogical(s)|ismember(s,[0 1])) % control saving
addParamValue(prs,'ifappend',true,@(s)islogical(s)|ismember(s,[0 1])) % if true, append new events to preexisting file
addParamValue(prs,'writing_behavior',[],@(s)isempty(s)|...
    ismember(s,{'append','replace','overwrite'})) % if not empty, controls writing behavior
addParamValue(prs,'FUNdefineEventsEpochs',@defineEventsEpochs_default)  % event/epoch definition file
addParamValue(prs,'filetype','behav',@(s)ismember(s,{'event' 'behav' 'stim'})) % 'event' ('behav') or 'stim'; determines what type of events to prealigning to
addParamValue(prs,'events',[],@(s)iscell(s)|isempty(s))  % events can be passed directly; definition file is not used
addParamValue(prs,'epochs',[],@(s)iscell(s)|isempty(s))  % epochs can be passed directly; definition file is not used
parse(prs,mycellids,varargin{:})
g = prs.Results;
if isempty(g.writing_behavior)  % keep backwards compatibility with 'ifappend' input argument
    if g.ifappend
        g.writing_behavior = 'append';
    else
        g.writing_behavior = 'overwrite';
    end
end
if ischar(mycellids)
    mycellids = {mycellids};    % convert to cell
end

% Get events and epochs
if isempty(g.events) && isempty(g.epochs)   % no events and epochs passed: use definition file
    [events,epochs] = feval(g.FUNdefineEventsEpochs);
else
    events = g.events;  % events and epochs were passed directly
    epochs = g.epochs;
end

% File name
switch g.filetype
    case {'behav','event'}
        FNAMESTR = 'EVENTSPIKES';
        EVFILE = 'TrialEvent';
    case 'stim'
        FNAMESTR = 'STIMSPIKES';
        EVFILE = 'StimEvent';
end

% Cell IDs
if strcmpi(mycellids,'all')  % if no cell IDs passed, run on all cells
    mycellids = listtag('cells');
end

% Load events
for i = 1:length(mycellids)
    cellid = mycellids(i);
    EV = loadcb(cellid,EVFILE);
    
    % Now we need to pass the Event file to extractEventSpikes
    [event_stimes,event_windows] = extractEventSpikes(cellid,events,EV);
    epoch_rates = extractEpochRates(event_stimes,event_windows,events,epochs);
    fname = cellid2fnames(cellid,FNAMESTR);
    
    % Save
    switch g.writing_behavior
        case 'append'  % try to append to existing file
            disp('Appending existing prealigned file.');
            x = load(fname);  % load existing alignment
            if strcmp(x.cellid,cellid)
                events = [x.events; events];  %#ok<AGROW>
                event_stimes = [x.event_stimes; event_stimes];  %#ok<NASGU,AGROW>
                event_windows = [x.event_windows; event_windows];  %#ok<NASGU,AGROW>
                epochs = [x.epochs; epochs];  %#ok<AGROW>
                epoch_rates = [x.epoch_rates; epoch_rates];  %#ok<NASGU,AGROW>
                if g.ifsave  % defensive strategy: save only if instructed to
                    save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates');
                end
            else
                disp('Mismatch in cellids. Append attempt is canceled. File not saved')
            end
        case 'replace'  % attempt to replace the redefined events
            disp('Appending existing prealigned file with replacing evnets.');
            x = load(fname);  % load existing alignment
            levents = events;   % 'local copies' of the newly computed variables - we need their names for the final ones
            levent_stimes = event_stimes;
            levent_windows = event_windows;
            lepochs = epochs;
            lepoch_rates = epoch_rates;
            events = x.events;   % start from the loaded variables and try to replace
            event_stimes = x.event_stimes;
            event_windows = x.event_windows;
            epochs = x.epochs;
            epoch_rates = x.epoch_rates;
            if strcmp(x.cellid,cellid)
                for iE = 1:size(levents,1)   % replace events
                    evinx = find(strcmp(events(:,1),levents{iE,1}));
                    if ~isempty(evinx)
                        events(evinx,:) = levents(iE,:);
                        event_stimes(evinx) = levent_stimes(iE);
                        event_windows(evinx) = levent_windows(iE);
                        if g.ifsave  % defensive strategy: save only if instructed to
                            save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates');
                        end
                    else
                        disp(['No matching event detected for ' events{iE} '. Replace attempt is canceled.'])
                    end
                end
                for iE = 1:size(lepochs,1)   % replace epochs
                    evinx = find(strcmp(epochs(:,1),lepochs{iE,1}));
                    if ~isempty(evinx)
                        epochs(evinx,:) = lepochs(iE,:);
                        epoch_rates(evinx) = lepoch_rates(iE);
                        if g.ifsave  % defensive strategy: save only if instructed to
                            save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates');
                        end
                    else
                        disp(['No matching epoch detected for ' epochs{iE} '. Replace attempt is canceled.'])
                    end
                end
            else
                disp('Mismatch in cellids. Append attempt is canceled. File not saved')
            end
        case 'overwrite'   % overwrite file
            if g.ifsave  % defensive strategy: save only if instructed to
                save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates');
            end
    end
end