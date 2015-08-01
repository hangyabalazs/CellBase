function cellbase_notes
%
% CellBase
%
% 2010-09-10
% TriggerName: Now the first field in the structure defineEventsEpochs will
% be the TriggerName and the second will be the TriggerEvent. This will
% enable proper plotting using variable events.

%% ADDANALYSIS
% SPR 2010-09-13
% added a NaN so that cellids where the analysis did not work will have
% NaNs instead of the entire analysis crashing.
%% DEFINEEVENTEPOCHS
% SPR 2010-09-13
% variable events
% for plotting rasters, the margin is [-2 0],and event names have suffix 2 
% for analysis of epoch rates, the margin is [0 0] and the event name does
% not have a suffix
% e.g. Use TriggerZoneHangOut2 to plot rasters and TriggerZoneHangOut to
% get the rate during the epoch.
% 2010-09-10
% *TRIGGERNAME* : Now the first field in the structure defineEventsEpochs will
% be the TriggerName and the second will be the TriggerEvent. This will
% enable proper plotting using variable events.
%% ADDNEWCELLS
% To add new cells to the database
% 1. Format the directory structure of the session to match everything else
% 2. run MakeTrialEvents2(sessionpath)
% 3. run MakeStimEvents2(sessionpath)
% 4. addnewcells('prealign')

% Since we do not keep any records of the epoch definitions we use, it is
% impossible to prealign spikes currently.
% This cn be acheived by opening up a cellid in the database and getting
% the event and epoch definition from that and running the prealignSpikes
% for the new cells using these definitions.
%% ADDCELL
% 1. this should run only after prealign. 
% 2. This also suffers from the fact that we may have changed the function since the time of running the analysis
% on all the cells in the database.
% (Need to store the .m file along with each property so we can have
% comparable analysis
% 3. Put in a try-catch so that if the property cannot be calculated, a NaN
% will be put in place.

%% ADDNEWSESSIONS
%% PREALIGNSPIKES
% Now you can add events and epochs as arguments
% e.g.
% prealignSpikes(allcells(1),'events',{behav_events},'epochs',{ebehav_epochs},'filetype','event')
% prealignSpikes(allcells(1),'events',{stim_events},'epochs',{stim_epochs},
% 'filetype','stim')
%% GET_MCLUST_WAVEFORMS4
% Now you can specify the spike times for which you want to calculate the
% spike waveforms as the second parameter
%% CALCULATE_CLUSEP
% script to calculate the ID and L-ratio and SNR for each cluster in the
% database.
% CAVEAT: It needs waveforms from all the channels otherwisr it will crash
% % CAVEAT removed, now it only considers channels with non-zero values
% The INPUT to function Create_CQ_File is the pwd, and it acts on all cells
% in the directory.
