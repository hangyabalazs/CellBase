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
mode = 'restrict';   % include only those sessions that contain the input cell IDs
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

% Plot parameters
gocolor = [0 0.8 0];   % color for go performace
nogocolor = [0.8 0 0];   % color for no-go performance
rtcolor = 'k';   % color for reaction time
ms = 8;   % marker size
lw = 2;   % line width

% Performance
[aSoundIntensity amGoPerformance amNoGoPerformance amGoReactionTime amNoGoReactionTime] = deal([]);
[aseGoPerformance aseNoGoPerformance aseGoReactionTime aseNoGoReactionTime] = deal([]);
for iM = 8:NumMice   % loop through mice
    animalID = mice{iM};   % current mouse
    sessions = findsession('animal',animalID);   % sessions of the current mouse
    NumSessions = size(sessions,1);   % number of sessions
    
    % Per animal
    [GoPerformance NoGoPerformance GoReactionTime NoGoReactionTime SoundIntensity] = deal([]);
    [ResponseRate_NoTone ResponseRate_NGTone ResponseRate_GoTone] = deal([]);
    for iS = 1:NumSessions   % loop through sessions
        sessionID = sessions{iS,2};   % current session
        cellIDs = findcell('rat',animalID,'session',sessionID);   % cells of the current session
        if isequal(mode,'restrict') && isempty(intersect(cellids,cellIDs))
            continue   % skip session if it does not contain any of the input cell IDs (e.g. no valid cluster from NB)
        end
        
        % Load trial events
        cellid = cellIDs{1};   % choose one cell from the session to load TrialEvents
        try
            TE = loadcb(cellid,'TrialEvents');
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
        
        % Session performance
        NumTrials = length(TE.Hit);
        
        % First restart
        longITIs = TE.ITIDistribution(1:end-1) > 1.8;   % long ITIs
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
        
        restartDist = restartitis(restartinx&longITIs);
        
        hit_inx = intersect(2:NumTrials-1,find(TE.Hit==1));
        fa_inx = intersect(2:NumTrials-1,find(TE.FalseAlarm==1));
        hit_inxm = unique([hit_inx-3 hit_inx-4 hit_inx-5]);
        hit_inxm(hit_inxm<1) = [];
        fa_inxm = unique([fa_inx-3 fa_inx-4 fa_inx-5]);
        fa_inxm(fa_inxm<1) = [];
        [GPhp NGPhp GORThp NGRThp] = psychplot_gonogo({animalID sessionID},'display',false,...
            'valid_trials',hit_inx+1);
        [GPhm NGPhm GORThm NGRThm] = psychplot_gonogo({animalID sessionID},'display',false,...
            'valid_trials',hit_inxm);
        [GPfp NGPfp GORTfp NGRTfp] = psychplot_gonogo({animalID sessionID},'display',false,...
            'valid_trials',fa_inx+1);
        [GPfm NGPfm GORTfm NGRTfm] = psychplot_gonogo({animalID sessionID},'display',false,...
            'valid_trials',fa_inxm);
        SI = unique(TE.SoundIntensity(1:end-1));   % different stimulus intesities
        if any(isnan(GPhp)) || any(isnan(GPhm)) || ...
            any(isnan(NGPhp)) || any(isnan(NGPhm)) % some conditions have no trials: not enough data
            continue
        end
        if ~isequal(length(SI),length(GPhp),length(GPhm)) % some conditions have no trials: not enough data
            continue
        end
        GoPerformance = [GoPerformance GPhp-GPhm]; %#ok<AGROW>
        NoGoPerformance =[NoGoPerformance NGPhp-NGPhm]; %#ok<AGROW>
        GoReactionTime =[GoReactionTime GORThp-GORThm]; %#ok<AGROW>
        NoGoReactionTime = [NoGoReactionTime NGRThp-NGRThm]; %#ok<AGROW>
        
        ResponseRate_NoTone = [ResponseRate_NoTone lickcount_notone]; %#ok<AGROW>
        ResponseRate_GoTone = [ResponseRate_GoTone lickcount_gtone]; %#ok<AGROW>
        ResponseRate_NGTone = [ResponseRate_NGTone lickcount_ngtone]; %#ok<AGROW>
        
        % Normalized sound intensity
        nSI = (SI - SI(1)) / (SI(end) - SI(1));
        SoundIntensity = [SoundIntensity nSI]; %#ok<AGROW>
    end
    
    mean(ResponseRate_NGTone)
    mean(ResponseRate_NoTone)
    mean(ResponseRate_GoTone)
    signrank(ResponseRate_NGTone,ResponseRate_NoTone)
    
    % Animal average
    if isempty(GoPerformance)  % no session was included in the analysis from this animal
        continue
    end
    [SI Ia Ib] = unique(SoundIntensity);   % different stimulus intesities
    NumInt = length(SI);   % number of different stimuli
    [mGoPerformance mNoGoPerformance mGoReactionTime mNoGoReactionTime] = ...
        deal(nan(1,NumInt));
    [seGoPerformance seNoGoPerformance seGoReactionTime seNoGoReactionTime] = ...
        deal(nan(1,NumInt));
    calc_average
    plot_average
    
    % Save
    if issave
        fnm = [resdir animalID '_PSY.fig'];
        saveas(H1,fnm)
        fnm = [resdir animalID '_PSY.jpg'];
        saveas(H1,fnm)
        fnm = [resdir animalID '_RT.fig'];
        saveas(H2,fnm)
        fnm = [resdir animalID '_RT.jpg'];
        saveas(H2,fnm)
        close([H1 H2])
    end
    
    NumInt = 3;   % pool intermediate intensities
    SoundIntensity2 = ones(size(SoundIntensity)) * 0.5;
    SoundIntensity2(SoundIntensity==0) = 0;
    SoundIntensity2(SoundIntensity==1) = 1;
    [SI Ia Ib] = unique(SoundIntensity2);   % different stimulus intesities
    [mGoPerformance mNoGoPerformance mGoReactionTime mNoGoReactionTime] = ...
        deal(nan(1,NumInt));
    [seGoPerformance seNoGoPerformance seGoReactionTime seNoGoReactionTime] = ...
        deal(nan(1,NumInt));
    calc_average
    plot_average
    
    % Save
    if issave
        fnm = [resdir animalID '_PSY2.fig'];
        saveas(H1,fnm)
        fnm = [resdir animalID '_PSY2.jpg'];
        saveas(H1,fnm)
        fnm = [resdir animalID '_RT2.fig'];
        saveas(H2,fnm)
        fnm = [resdir animalID '_RT2.jpg'];
        saveas(H2,fnm)
        close([H1 H2])
    end
    
    % Variables for grand average (pooled intermediate intensities)
%     amGoPerformance(end+1,:) = mGoPerformance; %#ok<AGROW>
%     amNoGoPerformance(end+1,:) = mNoGoPerformance; %#ok<AGROW>
%     amGoReactionTime(end+1,:) = mGoReactionTime; %#ok<AGROW>
%     amNoGoReactionTime(end+1,:) = mNoGoReactionTime; %#ok<AGROW>
%     aseGoPerformance(end+1,:) = seGoPerformance; %#ok<AGROW>
%     aseNoGoPerformance(end+1,:) = seNoGoPerformance; %#ok<AGROW>
%     aseGoReactionTime(end+1,:) = seGoReactionTime; %#ok<AGROW>
%     aseNoGoReactionTime(end+1,:) = seNoGoReactionTime; %#ok<AGROW>
    amGoPerformance = [amGoPerformance mGoPerformance]; %#ok<AGROW>
    amNoGoPerformance = [amNoGoPerformance mNoGoPerformance]; %#ok<AGROW>
    amGoReactionTime = [amGoReactionTime mGoReactionTime]; %#ok<AGROW>
    amNoGoReactionTime = [amNoGoReactionTime mNoGoReactionTime]; %#ok<AGROW>
    aseGoPerformance = [aseGoPerformance seGoPerformance]; %#ok<AGROW>
    aseNoGoPerformance = [aseNoGoPerformance seNoGoPerformance]; %#ok<AGROW>
    aseGoReactionTime = [aseGoReactionTime seGoReactionTime]; %#ok<AGROW>
    aseNoGoReactionTime = [aseNoGoReactionTime seNoGoReactionTime]; %#ok<AGROW>
    aSoundIntensity = [aSoundIntensity 0 0.5 1]; %#ok<AGROW>
end

% Grand average
NumInt = 3;   % intermediate intensities pooled
[SI Ia Ib] = unique(aSoundIntensity);   % different stimulus intesities
[mGoPerformance mNoGoPerformance mGoReactionTime mNoGoReactionTime] = ...
    deal(nan(1,NumInt));
[seGoPerformance seNoGoPerformance seGoReactionTime seNoGoReactionTime] = ...
    deal(nan(1,NumInt));
[GoPerformance NoGoPerformance GoReactionTime NoGoReactionTime] = ...
    deal(amGoPerformance,amNoGoPerformance,amGoReactionTime,amNoGoReactionTime);
calc_average
plot_average

% Save
if issave
    fnm = [resdir 'MEAN_PSY.fig'];
    saveas(H1,fnm)
    fnm = [resdir 'MEAN_PSY.jpg'];
    saveas(H1,fnm)
    fnm = [resdir 'MEAN_RT.fig'];
    saveas(H2,fnm)
    fnm = [resdir 'MEAN_RT.jpg'];
    saveas(H2,fnm)
end

% Overlay individual curves for each mouse
figure(H1)
hold on
plot(SI,reshape(amGoPerformance,3,length(amGoPerformance)/3),'Color',[0.7 1 0.7])
plot(SI,reshape(amNoGoPerformance,3,length(amNoGoPerformance)/3),'Color',[1 0.7 0.7])
figure(H2)
hold on
plot(SI,reshape(amGoReactionTime,3,length(amGoReactionTime)/3),'Color',[0.7 0.7 0.7])
if issave
    fnm = [resdir 'MEAN_PSY2.fig'];
    saveas(H1,fnm)
    fnm = [resdir 'MEAN_PSY2.jpg'];
    saveas(H1,fnm)
    fnm = [resdir 'MEAN_RT2.fig'];
    saveas(H2,fnm)
    fnm = [resdir 'MEAN_RT2.jpg'];
    saveas(H2,fnm)
    fnm = [resdir 'allvars.mat'];
    save(fnm)
end

GoTable = reshape(amGoPerformance,3,length(amGoPerformance)/3);
NoGoTable = reshape(amNoGoPerformance,3,length(amNoGoPerformance)/3);
for k = 1:3
    signrank(GoTable(k,:))
    signrank(NoGoTable(k,:))
end
if issave
    fnm = [resdir 'allvars.mat'];
    save(fnm)
end
keyboard

    % ---------------------------------------------------------------------
    function calc_average
        
        for idS = 1:NumInt
            inx = Ib == idS;
            mGoPerformance(idS) = mean(GoPerformance(inx));   % mean
            mNoGoPerformance(idS) = mean(NoGoPerformance(inx));
            mGoReactionTime(idS) = nanmean(GoReactionTime(inx));
            mNoGoReactionTime(idS) = nanmean(NoGoReactionTime(inx));
            seGoPerformance(idS) = nanse(GoPerformance(inx));   % SEM
            seNoGoPerformance(idS) = nanse(NoGoPerformance(inx));
            seGoReactionTime(idS) = nanse(GoReactionTime(inx));
            seNoGoReactionTime(idS) = nanse(NoGoReactionTime(inx));
        end
    end

    % ---------------------------------------------------------------------
    function plot_average
        
        % Plot animal average
        H1 = figure;   % performance
        set(gcf,'DefaultAxesTickDir','out')
        set(gcf,'DefaultAxesBox','off')
        xSI = linspace(0,1,length(SI));
        plot(xSI,mGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
            'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
        hold on
        errorbar(xSI,mGoPerformance,seGoPerformance,'Color',gocolor,'LineWidth',lw)
        plot(xSI,mNoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
            'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
        errorbar(xSI,mNoGoPerformance,seNoGoPerformance,'Color',nogocolor,'LineWidth',lw)
        ylim([-0.5 0.5])
        box off
        title('Performance')
        fstamp(animalID)
        
        H2 = figure;   % reaction time
        set(gcf,'DefaultAxesTickDir','out')
        set(gcf,'DefaultAxesBox','off')
        plot(xSI,mGoReactionTime,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
            'MarkerFaceColor',rtcolor,'Color',rtcolor,'MarkerSize',ms,'LineWidth',lw);
        errorbar(xSI,mGoReactionTime,seGoReactionTime,'Color',rtcolor,'LineWidth',lw)
        box off
        title('Reaction time')
        xlabel('Sound Pressure Level (dB)')
        fstamp(animalID)
    end
end