function [rate,dur,count] = extractEpochRates(event_stimes,event_windows,events,epochs)
%EXTRACTEPOCHRATES   Calculate spike rates in fixed or variable epochs.
%   [EPOCH_RATE,EPOCH_DUR,EPOCH_COUNT] = EXTRACTEPOCHRATES(EVENT_STIMES,EVENT_WINDOWS,EVENTS,EPOCHS)
%   EVENT_STIMES and EVENT_WINDOWS are 1xN cell arrays returned by 
%   EXTRACTEVENTSPIKES. EVENTS and EPOCHS are defined by DEFINEEVENTSEPOCHS.
%   It returns EPOCH_RATE, EPOCH_DUR, EPOCH_COUNT, each of which are 1xN 
%   cells arrays of 1xNumTrials arrays.
%
%   Using event windows, one can quickly calculate rates in different fixed
%   windows without re-extracting all the spike times. To use the whole of
%   a variable event (e.g. odor sampling or movement), use NaN for values
%   in that row.
%
%   See also EXTRACTEVENTSPIKES, DEFINEEVENTSEPOCHS and PREALIGNSPIKES.

%   Edit log: ZFM 10/04, AK 11/06, BH 6/27/11

% Initialize
NUMepochs = size(epochs,1);
rate = cell(1,NUMepochs);
dur = cell(1,NUMepochs);
count = cell(1,NUMepochs);

for iEP = 1:NUMepochs
    
    % Find the position of the event referenced in the epoch
    iEV = find(ismember(events(:,1),epochs{iEP,2}));
    
    % 'win0' and 'win1' are window times relative to trigger events      
    win0 = epochs{iEP,3}(1);
    win1 = epochs{iEP,3}(2);
    if isnan(win0) || isnan(win1)
        VARWINDOW = 1;
    else
        VARWINDOW = 0;
    end
    
    for iT = 1:length(event_stimes{iEV})
        
        % Window for variable length epochs taken per trial
        if VARWINDOW
            w0 = event_windows{iEV}(1,iT);
            w1 = event_windows{iEV}(2,iT);
        else
            w0 = win0;
            w1 = win1;
        end
        
        % Check if w0 or w1 is outside the windows used to extract the spike times
        if isnan(w0) || isnan(w1) || w0 < event_windows{iEV}(1,iT) || w1 > event_windows{iEV}(2,iT)
            count{iEP}(iT) = NaN;
            dur{iEP}(iT) = NaN;
            rate{iEP}(iT) = NaN;
        else
            count{iEP}(iT) = length(find(event_stimes{iEV}{iT} >= w0 & event_stimes{iEV}{iT} < w1));
            dur{iEP}(iT) = w1 - w0;
            rate{iEP}(iT) = count{iEP}(iT) / dur{iEP}(iT);
        end
    end   % iT
end   % iEP