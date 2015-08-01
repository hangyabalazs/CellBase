function psych_gonogo(cellids)
%PSYCH_GONOGO   Average performance and reaction time.
%   PSYCH_GONOGO(CELLIDS) calculates average performance and reaction time
%   for all sessions corresponding to CELLIDS, per animal. Sessions that
%   contain light-stimulation are excluded.
%
%   See also PSYCHPLOT_GONOGO.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Sept-2013

%   Edit log: BH 9/27/13

% Input arguments
mode = 'restrict';   % include only those sessions that contain the input cell IDs

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
for iM = 1:NumMice   % loop through mice
    animalID = mice{iM};   % current mouse
    sessions = findsession('animal',animalID);   % sessions of the current mouse
    NumSessions = size(sessions,1);   % number of sessions
    
    % Per animal
    [GoPerformance NoGoPerformance GoReactionTime NoGoReactionTime SoundIntensity] = deal([]);
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
        [GP NGP GORT NGRT] = psychplot_gonogo({animalID sessionID},'display',false);
        SI = unique(TE.SoundIntensity(1:end-1));   % different stimulus intesities
        GoPerformance = [GoPerformance GP]; %#ok<AGROW>
        NoGoPerformance =[NoGoPerformance NGP]; %#ok<AGROW>
        GoReactionTime =[GoReactionTime GORT]; %#ok<AGROW>
        NoGoReactionTime = [NoGoReactionTime NGRT]; %#ok<AGROW>
        
        % Normalized sound intensity
        nSI = (SI - SI(1)) / (SI(end) - SI(1));
        SoundIntensity = [SoundIntensity nSI]; %#ok<AGROW>
    end
    
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
    for iS = 1:NumInt
        inx = Ib == iS;
        mGoPerformance(iS) = mean(GoPerformance(inx));   % mean
        mNoGoPerformance(iS) = mean(NoGoPerformance(inx));
        mGoReactionTime(iS) = mean(GoReactionTime(inx));
        mNoGoReactionTime(iS) = mean(NoGoReactionTime(inx));
        seGoPerformance(iS) = std(GoPerformance(Ib==iS)) / sqrt(sum(inx));   % SEM
        seNoGoPerformance(iS) = std(NoGoPerformance(Ib==iS)) / sqrt(sum(inx));
        seGoReactionTime(iS) = std(GoReactionTime(Ib==iS)) / sqrt(sum(inx));
        seNoGoReactionTime(iS) = std(NoGoReactionTime(Ib==iS)) / sqrt(sum(inx));
    end
    
    % Plot animal average
    figure   % performance
    set(gcf,'DefaultAxesTickDir','out')
    set(gcf,'DefaultAxesBox','off')
    xSI = linspace(0,1,length(SI));
    plot(xSI,mGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
    errorbar(xSI,mGoPerformance,seGoPerformance,'Color',gocolor,'LineWidth',lw)
    hold on
    plot(xSI,mNoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
        'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
    errorbar(xSI,mNoGoPerformance,seNoGoPerformance,'Color',nogocolor,'LineWidth',lw)
    ylim([0 1])
    box off
    title('Performance')
    fstamp(animalID)
    
    figure   % reaction time
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