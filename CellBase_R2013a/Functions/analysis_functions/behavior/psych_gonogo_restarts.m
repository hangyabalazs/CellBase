function psych_gonogo_restarts(cellids,issave)
%PSYCH_GONOGO_UPDATING   Average performance and reaction time.
%   PSYCH_GONOGO_UPDATING(CELLIDS,ISSAVE) calculates average performance
%   and reaction time for all sessions corresponding to CELLIDS, per
%   animal. Sessions that contain light-stimulation are excluded. Trials
%   before and after hits or false alarms are compared to check for
%   updating effects.
%   Input parameters:
%       CELLIDS - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, all cells are
%           selected from the CellBase
%       ISSAVE - controls saving
%   Plots:
%       1. Average performance and reaction time per mouse for the sessions
%       in which the cells (CELLIDS) were recoreded. The sound intensities
%       are normalized between 0 and 1 within each session. After averaging
%       according to normalized intensities, the intensities are evenly
%       distributed for plotting.
%       2. In another within-mouse average, all intermediate intensities
%       are pooled within a session. For plotting, these intensities are
%       displayed in the middle.
%       3. Grand average: average across mice from the second type of
%       within-mouse averages.
%
%   See also PSYCHPLOT_GONOGO.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Sept-2013

%   Edit log: BH 9/27/13, 11/1/13

% Pass the control to the user in case of error
dbstop if error

% Input arguments
mode = 'nonrestrict';   % include only those sessions that contain the input cell IDs
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1 || isempty(cellids)
    loadcb   % load CellBase
    cellids = CELLIDLIST;
else
    if isnumeric(cellids)
        loadcb   % load CellBase
        cellids = CELLIDLIST(cellids);   % index set to CELLIDLIST
    elseif ischar(cellids)
        cellids = {cellids};   % only one cellID passed
    elseif iscellstr(cellids)
        % list of cell IDs
    else
        error('psych_gonogo:inputArg','Unsupported format for cell IDs.')
    end
end

% Directories
global DATAPATH
resdir = fullfile(DATAPATH,'NB','average_performance',filesep,'restarts',filesep);

% Animals
mice = listtag('animal');
NumMice = length(mice);   % number of animals
basedir = getpref('cellbase','datapath');

% Performance
for iM = 11:NumMice   % loop through mice
    animalID = mice{iM};   % current mouse
    sessions = findallsession('animal',animalID);   % sessions of the current mouse
%     sessions = sessions(size(sessions,1)-10:end,:);
    NumSessions = size(sessions,1);   % number of sessions
    
    % Per animal
    [ResponseRate_NoTone ResponseRate_NGTone ResponseRate_GoTone] = deal([]);
    for iS = 1:NumSessions   % loop through sessions
        sessionID = sessions{iS,2};   % current session
        cellIDs = findcell('rat',animalID,'session',sessionID);   % cells of the current session
        if isequal(mode,'restrict') && isempty(intersect(cellids,cellIDs))
            continue   % skip session if it does not contain any of the input cell IDs (e.g. no valid cluster from NB)
        end
        
        % Load trial events
        datapath = fullfile(basedir,animalID,sessionID);
        try
            TE = load([datapath filesep 'TE.mat']);
        catch %#ok<CTCH>
            disp([sessionID ': No TrialEvents file.'])
            continue
        end
        
        % Exclude sessions with light stimulation
        if isfield(TE,'LightStimulation2')   % exclude sessions with light-stimulation (affects only one session of NB CellBase)
            lighton_trials = find(TE.LightStimulation2==1,1);
            if ~isempty(lighton_trials)
                disp([animalID ' ' sessionID ': behavior session with light-stimulation - excluded.'])
                continue
            end
        end
        
        if ~isequal(nanmin(TE.StimulusDuration),20)   % exclude early sessions with no hard trials
            continue
        end
        
        % Session performance
        NumTrials = length(TE.Hit);
        
        % First restart
        longITIs = TE.ITIDistribution(1:end-1) > 1.8;   % long ITIs
%         shortITIs = TE.ITIDistribution(1:end-1) < 1.2;
        shortITIs = TE.ITIDistribution(1:end-1) < 1.2 & ...
            TE.StimulusDuration(1:end-1)==nanmax(TE.StimulusDuration);   % short ITIs, high sound intensity
        
        restartinx = cellfun(@(s)length(s)>1,TE.ITIBegins(1:end-1));   % more than one ITIs
        restartbegins = cellfun(@(s)s(1),TE.ITIBegins(1:end-1));
        restartends = cellfun(@(s)s(1),TE.ITIEnds(1:end-1));
        restartitis = restartends - restartbegins;   % length of first ITI in trial (full ITI if no restart)
        
        nogos = TE.CorrectRejection(1:end-1) == 1 | TE.FalseAlarm(1:end-1) == 1;   % No-go trials
        fainx = TE.FalseAlarm(1:end-1) == 1;   % False Alarms
        nonnaninx = nogos & shortITIs & (restartinx | fainx);   % lick to (loud) no-go in short ITI trials
        allinx = nogos & shortITIs;
        temprestartitis = restartitis;
        temprestartitis(allinx&~nonnaninx) = NaN;   % nan out non-restarted correct rejections
        nogo_restarts = temprestartitis(allinx);
        lickcount_ngtone = sum(nogo_restarts>TE.ITIDistribution(allinx)&...
            nogo_restarts<TE.ITIDistribution(allinx)+0.6) / length(nogo_restarts);   % short ITI false alarm rate
        
        flongitis = find(longITIs);   % long ITIs
        randsamp = randi(sum(longITIs),[1 length(nogo_restarts)]);   % random sample matched to no-go licks
        long_restarts = restartitis(flongitis(randsamp));
        long_restarts(long_restarts>TE.ITIDistribution(allinx)+0.6) = NaN;
        lickcount_notone = sum(long_restarts>TE.ITIDistribution(allinx)&...
            long_restarts<TE.ITIDistribution(allinx)+0.6) / length(long_restarts);   % count restars in the no-go windows
        
        gos = TE.Hit(1:end-1) == 1 | TE.Miss(1:end-1) == 1;
        hitinx = TE.Hit(1:end-1) == 1;
        nonnaninx = gos & shortITIs & (restartinx | hitinx);
        allinx = gos & shortITIs;
        temprestartitis = restartitis;
        temprestartitis(allinx&~nonnaninx) = NaN;
        go_restarts = temprestartitis(allinx);
        lickcount_gtone = sum(go_restarts>TE.ITIDistribution(allinx)&...
            go_restarts<TE.ITIDistribution(allinx)+0.6) / length(go_restarts);
        
        ResponseRate_NoTone = [ResponseRate_NoTone lickcount_notone]; %#ok<AGROW>
        ResponseRate_GoTone = [ResponseRate_GoTone lickcount_gtone]; %#ok<AGROW>
        ResponseRate_NGTone = [ResponseRate_NGTone lickcount_ngtone]; %#ok<AGROW>
    end
    
    mean(ResponseRate_NGTone)
    mean(ResponseRate_NoTone)
    mean(ResponseRate_GoTone)
    signrank(ResponseRate_NGTone,ResponseRate_NoTone)
    1;
end 