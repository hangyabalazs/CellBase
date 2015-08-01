function NUM_ADDED = addnewcells(varargin)
%ADDNEWCELLS   Add new cells to CellBase.
%   NM = ADDNEWCELLS adds new cells found in Cellbase directory structure
%   to CellBase calling FINDALLCELLS. It calls ADDCELL, which also performs
%   all the analyses. The number of newly added cells is returned.
%
%   NM = ADDNEWCELLS('PREALIGN',TRUE) calls PREALIGNSPIKES.
%
%   NM = ADDNEWCELLS('QUIET',TRUE) suppresses message displays.
%
%   NM = ADDNEWCELLS('DIR',DR) looks for new cells in DR. DR should contain
%   the path downstream of the CellBase directory, e.g.
%   addnewcells('dir','n046\130101a\')
%
%   See also FINDALLCELLS and ADDCELL.

%   Edit log: AK 10/06, BH 3/21/11, 4/30/14

% Default arguments
prs = inputParser;
addParamValue(prs,'prealign',false,@(s)islogical(s)|ismember(s,[0 1]))   % control whether to run prealignSpikes
addParamValue(prs,'quiet',false,@(s)islogical(s)|ismember(s,[0 1]))   % control command line messages
addParamValue(prs,'force',false,@(s)islogical(s)|ismember(s,[0 1]))   % force mode disabled for now
addParamValue(prs,'dir','',@ischar)   % restrict addnewcells to a directory
parse(prs,varargin{:})
g = prs.Results;

% Get CellBase preferences
clear global CELLIDLIST ANALYSES TheMatrix   % refresh the globals
cellbase_datapath = getpref('cellbase','datapath');

% Find all cells and sort the new ones
if isempty(g.dir)
    all_cellids = findallcells;
else
    all_cellids = findallcells(g.dir);
end
old_cellids = listtag('cells');
clear global CELLIDLIST ANALYSES TheMatrix   % refresh the globals

if isempty(old_cellids)
    new_cellids = all_cellids;  % setdiff fails for empty set
else
    new_cellids = setdiff(all_cellids,old_cellids);
end
NUM_newcellids = length(new_cellids);
global NEWCELLIDS
NEWCELLIDS = new_cellids;   % enables adding analysis (with addanalysis) that runs iscellid

% Add cells; prealign if called that way
NUM_ADDED = 0;
if NUM_newcellids > 0
    if ~g.quiet
        disp(sprintf('%d new cells found.',NUM_newcellids))
    end
    
    if g.prealign
        for iOld = 1:length(old_cellids), % get the latest events and epochs in the database
            try
                EVENTSPIKES = loadcb(old_cellids(iOld),'EventSpikes');
            catch
                EVENTSPIKES = [];
            end
            if ~isempty(EVENTSPIKES)
                behav_events = EVENTSPIKES.events;
                behav_epochs = EVENTSPIKES.epochs;
                break
            else
            end
        end
        for iOld = 1:length(old_cellids)  % get the latest events and epochs in the database
            try
                STIMSPIKES = loadcb(old_cellids(iOld),'StimSpikes');
            catch
                STIMSPIKES = [];
            end
            if ~isempty(STIMSPIKES)
                stim_events = STIMSPIKES.events;
                stim_epochs = STIMSPIKES.epochs;
                break
            else
            end
        end
        for iNuCell = 1:length(new_cellids)
            try
                prealignSpikes(new_cellids(iNuCell),'events',{behav_events},'epochs',{behav_epochs},'filetype','behav')
            catch
            end
            try
                prealignSpikes(new_cellids(iNuCell),'events',{stim_events},'epochs',{stim_epochs},'filetype','stim')
            catch
            end
        end
    end
    for iC = 1:NUM_newcellids
        OK = addcell(new_cellids{iC},'quiet');  % add cells
        if OK
            NUM_ADDED = NUM_ADDED + 1;
        end
    end   %iC
else   %NUM_newcellids
    if ~g.quiet
        disp('No new cells found.')
    end
end