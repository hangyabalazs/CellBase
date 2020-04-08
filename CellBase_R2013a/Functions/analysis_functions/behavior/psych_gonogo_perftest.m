function psych_gonogo_learningcurve(cellids,issave)
%PSYCH_GONOGO_LEARNINGCURVE   Average learning curve.
%   PSYCH_GONOGO_LEARNINGCURVE(CELLIDS,ISSAVE) calculates average learning
%   curve corresponding to CELLIDS, per animal. The easiest stimulus
%   condition is used. Trials after the animals' last response are
%   excluded.
%   Input parameters: 
%       CELLIDS - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, all cells are
%           selected from the CellBase
%       ISSAVE - controls saving
%   Plots:
%       1. Average learning curve per mouse for the animals in which the
%       cells (CELLIDS) were recoreded.
%       2. Grand average: average across mice.
%
%   See also PSYCHPLOT_GONOGO.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   3-Nov-2013

%   Edit log: 11/3/13

%   Note for NB mice: Animal 15 is excluded due to missing data and animal
%   40 is excluded due to early stimulation sessions. 

% Pass the control to the user in case of error
dbstop if error

% Input arguments
mode = 'restrict';   % include only those sessions that contain the input cell IDs
error(nargchk(0,2,nargin))
if nargin < 2
    issave = false;
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
        error('psych_gonogo_learningcurve:inputArg','Unsupported format for cell IDs.')
    end
end

% Directories
global DATAPATH
resdir = fullfile(DATAPATH,'NB','perftest',filesep);

% Animals
mice = listtag('animal');
mice(strcmp(mice,'n015')|strcmp(mice,'n040')) = [];   % exclude due to stim. (n040) or missing behav. data (n015)
NumMice = length(mice);   % number of animals

% Plot parameters
gocolor = [0 0.8 0];   % color for go performace
nogocolor = [0.8 0 0];   % color for no-go performance
rtcolor = 'k';   % color for reaction time
ms = 8;   % marker size
lw = 2;   % line width

% Performance
[aGoPerformance aNoGoPerformance Ps] = deal(nan(50,200));
figure
hold on
for iM = 1:NumMice   % loop through mice
    animalID = mice{iM};   % current mouse
%     sessions = findsession('animal',animalID);   % sessions of the current mouse
    inpdir = fullfile(DATAPATH,'\NB\_behavior\',[animalID(1) 'b' animalID(2:end)]);
    psi = dir(inpdir);
    sessionIDs = {psi(3:end).name};   % all sessions, not only those having cells
    sessionIDs = cellfun(@(s)s(1:7),sessionIDs,'UniformOutput',false);
    sessionIDs = unique(sessionIDs)';
    NumSessions = size(sessionIDs,1);   % number of sessions
    
    cellIDs = findcell('rat',animalID);   % cells of the current mouse
    if isequal(mode,'restrict') && isempty(intersect(cellids,cellIDs))
        continue   % skip mouse if it does not contain any of the input cell IDs (e.g. no valid cluster from NB)
    end
    
    % Per animal
    [GoPerformance NoGoPerformance] = deal([]);
    P = [];
    for iS = 1:NumSessions   % loop through sessions
        sessionID = sessionIDs{iS};   % current session
                
        % Load trial events
        sessiondr = fullfile(getpref('cellbase','datapath'),animalID,sessionID,filesep);
        TE = load(fullfile(sessiondr,'TE.mat'));
        
        NumTrials = length(TE.TrialStart);
        wn = 20;   % window size (number of trials)
        [hits,fas,gos,nogos,hitrate,farate,restartrate] = deal(nan(NumTrials-1,1));
        for iT = 1:NumTrials-1     % last trial may be incomplete
            
            % Outcome
            if iT >= wn
                hits(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT)));   % number of hits
                fas(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT)));   % number of false alarms
                gos(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT))|~isnan(TE.Miss(iT-wn+1:iT)));   % number of go tones
                nogos(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT))|~isnan(TE.CorrectRejection(iT-wn+1:iT)));   % number of no-go tones
                
                %             restartrate(iT) = sum(numrestarts(iT-wn+1:iT));   % number of restarts
            end
            
            hitrate(iT) = hits(iT) / gos(iT);   % hit rate
            farate(iT) = fas(iT) / nogos(iT);   % false alarm rate
            
        end
        
        startinx = find(hitrate>=0.7,1,'first');   % engagement: from first hit
        endinx = find(hitrate>=0.7,1,'last');
        if isempty(startinx) || isempty(endinx)
            p = NaN;
            P = [P p];
            continue
        end
        ln(iS) = endinx - startinx;   % period of engagement
        
        ehits = sum(~isnan(TE.Hit(startinx:endinx)));   % number of hits in the engaged period
        efas = sum(~isnan(TE.FalseAlarm(startinx:endinx)));   % number of false alarms in the engaged period
        egos = sum(~isnan(TE.Hit(startinx:endinx))|~isnan(TE.Miss(startinx:endinx)));   % number of go tones in the engaged period
        enogos = sum(~isnan(TE.FalseAlarm(startinx:endinx))|~isnan(TE.CorrectRejection(startinx:endinx)));   % number of no-go tones in the engaged period
        ehitrate(iS) = ehits / egos;   % hit rate
        efarate(iS) = efas / enogos;   % false alarm rate
        
%         figure
%         plot(hitrate,'g')
%         hold on
%         plot(farate,'r')
        
        % Chi-square test
        [h p] = b_chi2test2([ehits egos-ehits],[efas enogos-efas],0.05)   % chi square test
        P = [P p];   % test for sign. perf.
        
        % Session performance
%         [GP NGP] = psychplot_gonogo({animalID sessionID},...
%             'lastresponsestop',true,'display',false);
%         GoPerformance = [GoPerformance GP(end)]; %#ok<AGROW>
%         NoGoPerformance =[NoGoPerformance NGP(end)]; %#ok<AGROW>
    end
    nm = length(P);
    Ps(iM,1:nm) = P;
    plot(P)
%     aGoPerformance(iM,1:nm) = GoPerformance;
%     aNoGoPerformance(iM,1:nm) = NoGoPerformance;
%     
%     % Plot
%     H1 = figure;
%     set(gcf,'DefaultAxesTickDir','out')
%     set(gcf,'DefaultAxesBox','off')
%     plot(GoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
%         'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
%     hold on
%     plot(NoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
%         'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
%     ylim([0 1])
%     box off
%     title('Learning curve')
%     xlabel('%lick')
%     ylabel('Session')
%     fstamp(animalID)
%         
%     % Save
%     if issave
%         fnm = [resdir animalID '_LearningCurve.fig'];
%         saveas(H1,fnm)
%         fnm = [resdir animalID '_LearningCurve.jpg'];
%         saveas(H1,fnm)
%         close(H1)
%     end
end
keyboard


% % Grand average
% mGoPerformance = nanmean(aGoPerformance);   % mean
% mNoGoPerformance = nanmean(aNoGoPerformance);
% seGoPerformance = nanse(aGoPerformance);   % SEM
% seNoGoPerformance = nanse(aNoGoPerformance);
% 
% % Plot
% H1 = figure;   % performance
% set(gcf,'DefaultAxesTickDir','out')
% set(gcf,'DefaultAxesBox','off')
% plot(mGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
%     'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
% hold on
% E = errorbar(mGoPerformance,seGoPerformance,'Color',gocolor,'LineWidth',lw);
% errorbar_tick(E,0)   % eliminate horizontal line from errorbar
% plot(mNoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
%     'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
% E = errorbar(mNoGoPerformance,seNoGoPerformance,'Color',nogocolor,'LineWidth',lw);
% errorbar_tick(E,0)   % eliminate horizontal line from errorbar
% ylim([0 1])
% box off
% title('Learning curve')
% xlabel('%lick')
% ylabel('Session')
%     
% % Save
% if issave
%     fnm = [resdir 'MEAN_LC.fig'];
%     saveas(H1,fnm)
%     fnm = [resdir 'MEAN_LC.jpg'];
%     saveas(H1,fnm)
%     fnm = [resdir 'MEAN_LC.mat'];
%     save(fnm,'aGoPerformance','aNoGoPerformance')
% end
% 
% % Overlay individual curves for each mouse
% figure(H1)
% hold on
% plot(aGoPerformance','Color',[0.7 1 0.7])
% plot(aNoGoPerformance','Color',[1 0.7 0.7])
% if issave
%     fnm = [resdir 'MEAN_LC2.fig'];
%     saveas(H1,fnm)
%     fnm = [resdir 'MEAN_LC2.jpg'];
%     saveas(H1,fnm)
% end

% -------------------------------------------------------------------------
function [h,p,stat] = b_chi2test2(n_d1,n_d2,alpha)
% for binned data

if nargin < 3
    alpha = 0.05;
end

ld1 = sum(n_d1);
ld2 = sum(n_d2);

stat = ld1 * ld2 * sum((n_d1/ld1-n_d2/ld2).^2./(n_d1+n_d2));
q = chi2cdf(stat,length(n_d1)-1);
if q > (1 - alpha)
    h = 1;
else
    h = 0;
end
p = 1 - q;