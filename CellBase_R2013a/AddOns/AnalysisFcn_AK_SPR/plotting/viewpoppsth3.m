function handles = viewpoppsth3(cellids,partition,varargin)
% Allows plotting dashed condition for all psths based on some criterion event.
%   VIEWPOPPSTH3
%
% SPR 2012-02-02

if (nargin < 1)
    help viewpoppsth3
    return
end

% check if cellid is valid  validcellid(cellid,{'list'}) ~= 1

if ~isempty(varargin)
    try, g = struct(varargin{:});
    catch, error('Argument error in the {''param'', value} sequence'); end;
end;

% assign defaults if no arguments passed
try, g.TriggerName;    catch, g.TriggerName  = 'WaterPokeIn'; end;
try, g.OdorPairID;      catch, g.OdorPairID    = 1;  end;
try, g.Normalization;   catch, g.Normalization = 'max'; end;
try, g.window;          catch, g.window        = [-0.5 1];  end;
try, g.dt;              catch, g.dt            = 0.01; end;
try, g.sigma;           catch, g.sigma         = 0.02; end;
try, g.plot;            catch, g.plot          = 'on'; end;
try, g.FigureNum;       catch, g.FigureNum     = 1; end;
try, g.ValidTrials;     catch, g.ValidTrials   = ''; end;
try, g.Overlay;         catch, g.Overlay       = 'off'; end;
try, g.BurstPSTH;       catch, g.BurstPSTH     = 'off'; end;
try, g.PSTHstd;       catch, g.PSTHstd     = 'off'; end;
try, g.PlotDashedEvent;       catch, g.PlotDashedEvent     = ''; end;
try, g.PlotDashedCondition;       catch, g.PlotDashedCondition     = 'min'; end;
try, g.DashedLineStyle;             catch, g.DashedLineStyle='--'; end;
try, g.PrevEvent;  catch, g.PrevEvent = ''; end
try, g.NextEvent;  catch, g.NextEvent = ''; end
try, g.spikerestrict;  catch, g.spikerestrict= 'none'; end


% test arguments for consistency
switch lower(g.plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;
if ~isempty(g.ValidTrials)
    g.ValidTrials = strcat('& ',g.ValidTrials);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumCells = length(cellids);
% NORMstimrateC = zeros(NumCells,8)*NaN;
% NORMstimrateE = zeros(NumCells,8)*NaN;
margin = g.sigma*3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array

% TE = loadcb(cellids{1},'Events');
TE = loadcb(cellids{1},'Trialevent');
SP = loadcb(cellids{1},'EventSpikes');
TEfields = fieldnames(TE);

Npart = sum(ismember(lower(TEfields),lower(partition)));
if Npart == 0,
    Npart = size(partition,2);
end
if ismember('mixturediff',lower(partition))
    Npart=Npart+3;
end
if ismember('odorratio',lower(partition))
    Npart=Npart+6;
end

NPSTH=zeros(NumCells,Npart,length(time))*NaN;
g.PlotDashedTime = nan(NumCells,Npart);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d cells to process.\n',NumCells)
for iCELL = 1:NumCells
    
    print_progress(iCELL,round(NumCells/100),5);
    
    cellid = cellids{iCELL};
    % TE = loadcb(cellid,'Events');
    TE = loadcb(cellid,'Trialevent');
    ST = loadcb(cellid,'EVENTSPIKES');
    % ST = loadcb(cellid,'STIMES');
    
    trigger_pos = findcellstr(ST.events(:,1),g.TriggerName);
    
    if (trigger_pos == 0)
        error('Trigger variable not found');
    else
         g.TriggerEvent=SP.events{trigger_pos,2};
    end
    
    alltrials = 1:size(ST.event_stimes{1},2);
    
    stimes  = ST.event_stimes{trigger_pos}(alltrials);
    
    % get windows for previous event.
    
        % restrict spikes from flanking events.
    switch g.spikerestrict
        case 'none'
            windows = repmat([g.window(1)-margin g.window(2)+margin],size(stimes,2),1);
        case 'both'
            windows = [trialevents2relativetime(TE,g.TriggerEvent,{g.PrevEvent})'...
                trialevents2relativetime(TE,g.TriggerEvent,{g.NextEvent})'];
        case 'pre'
            windows = [trialevents2relativetime(TE,g.TriggerEvent,{g.PrevEvent})'...
                (g.window(2)+margin)*(ones(size(stimes,2),1))];
        case 'post'
              windows = [(g.window(1)-margin)*(ones(size(stimes,2),1)) ...
                  trialevents2relativetime(TE,g.TriggerEvent,{g.NextEvent})'];
    end
%     margins = [0 0];

    
%     
%     if ~isempty(g.VariableEvent),
%         winevent_pos = findcellstr(ST.events(:,1),g.VariableEvent);
%     end
%     if ~isempty(winevent_pos),
%         windows = ST.event_windows{winevent_pos}(:,alltrials);
%     else
%         windows = ST.event_windows{trigger_pos}(:,alltrials);
%     end
%     
%     windows = -flipud(windows);
%     windows(2,:) = g.window(2);
%    
    if strcmp(g.BurstPSTH,'on'),
        stimes = detect_bursts(stimes);
    end
    
    %%% MAKE THE MAIN RASTER
%     binraster = stimes2binraster(stimes,time,g.dt);
    binraster = stimes2binraster(stimes,time,g.dt,windows);
    
    %compare = {'Correct','Error','ChoiceLeft','ChoiceRight','LeftCorrect','LeftError','RightCorrect','RightError','Stimulus'};
    
    
    clear COMP;
    iP = 1; iC = 1;
    while length(partition) >= iP
        
        switch lower(partition{iP});
            case 'odorratio'
                Ratios = setdiff(unique(TE.OdorRatio),50);
                NumRatio = length(Ratios);
                for i=1:NumRatio
                    COMP(iC+i-1,:) = eval(['TE.OdorRatio==' num2str(Ratios(i))]);
                    Plegend{iC+i-1} = num2str(Ratios(i));
                end
                iC=iC+i+1;
            case  'mixturediff'
                Ratios = setdiff(unique(TE.OdorRatio),50);
                NumRatio = length(Ratios);
                H = NumRatio / 2;
                COMP(iC,:) = (TE.OdorRatio == Ratios(H)) | (TE.OdorRatio == Ratios(H+1));
                Plegend{iC} = 'Diff';
                iC=iC+1;
                if NumRatio > 4
                    COMP(iC,:) = (TE.OdorRatio == Ratios(H-1)) | (TE.OdorRatio == Ratios(H+2));
                    COMP(iC+1,:) = (TE.OdorRatio == Ratios(H-2)) | (TE.OdorRatio == Ratios(H+3));
                    Plegend{iC} = 'Med';
                    Plegend{iC+1} = 'Easy';
                    iC=iC+2;
                else
                    COMP(iC+1,:) = (TE.OdorRatio == Ratios(H-1)) | (TE.OdorRatio == Ratios(H+2));
                    Plegend{iC+1} = 'Easy';
                    iC=iC+2;
                end
            case 'all'
                COMP(iC,:) = ones(1,length(TE.TrialStart));
                Plegend{iC} = partition{iP};
                iC=iC+1;
            case lower(TEfields)
                COMP(iC,:) = eval(['TE.' partition{iP}]);
                Plegend{iC} = partition{iP};
                iC=iC+1;
        end
        iP=iP+1;
        if ~exist('COMP','var'),
                [COMP, Plegend] = partition_trials(TE,partition);
                continue
        end
    end %while partition
    
    
    % if g.OdorPairID == 0
    %     valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.ValidTrials));
    % else
    %     valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.OdorPairID,g.ValidTrials));
    % end
    
    valid_trials = 1:length(TE.TrialStart);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(find(valid_trials)) > 10
        [psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMP,valid_trials);
    else
        spsth    = NaN;
        spsth_se = NaN;
    end
    numpsth = size(spsth,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(ismember(partition,'OdorRatio'))
        if NumRatio <= 4
            spsth(3,:) = spsth(2,:)*NaN;
        end
    end
    
    NORMbaseline = 0;
    switch lower(g.Normalization)
        case 'none'
            NORMfactor = 1;
        case 'mean'
            NORMfactor =  nanmean(spsth(:));
        case 'median'
            NORMfactor =  nanmedian(spsth(:));
        case 'max'
            NORMfactor =  nanmax(spsth(:));
        case 'min'
            NORMbaseline = nanmean(spsth(:));
            NORMfactor =  -nanmin(spsth(:)-NORMbaseline);
        case 'zscore'
            NORMfactor =  nanstd(spsth(:));
            NORMbaseline = nanmean(spsth(:));
            %     case 'maxrate'
            %         NORMfactor = max(EpochRate(union(trialsC,trialsE)));
            %     case 'maxratecorr'
            %         NORMfactor = max(EpochRate(trialsC));
            %     case 'maxrateerr'
            %         NORMfactor = max(EpochRate(trialsE));
            %     case 'perc95'
            %          NORMfactor = prctile(EpochRate(trialsE),90);
            %     case 'meanrate'
            %         NORMfactor = nanmean(EpochRate(union(trialsC,trialsE)))
        otherwise
            sprintf('viewpoppsth: Normalization method not found. NORMfactor = 1; \n');
            NORMfactor = 1;
    end
    
    
    NPSTH(iCELL,1:size(spsth,1),1:size(spsth,2)) = (spsth-NORMbaseline)/ NORMfactor;
    
    % Plot Dashed Condition on POPPSTH
    if ~isempty(g.PlotDashedEvent),
        %junk = valid_trials;
        % i.e. TrialStart is an actual timestamp
        if nansum(TE.TrialStart)~=0,
            DashedTimes = eval(['TE.' g.PlotDashedEvent]);
            % if TrialStart is a dummy variable set to 0
        else
            % if DashedEvent is before TriggerEvent
            DashedTimes = trialevents2relativetime(TE,g.TriggerEvent,{g.PlotDashedEvent});
% if nanmean(TE.(g.TriggerEvent)  - TE.(g.PlotDashedEvent) ) > 0,
%                 DashedTimes = TE.(g.PlotDashedEvent) - TE.(g.TriggerEvent);
%                 % if DashedEvent is after TriggerEvent
%             else
%                 DashedTimes = TE.(g.TriggerEvent) - TE.(g.PlotDashedEvent);
%             end
            %         DashedTimes=SP.event_windows{trigger_pos}(2,valid_trials);
        end
%         DashedTimes = DashedTimes(valid_trials);
        
        % get PlotDashedTime for each psth
        for iD = 1:numpsth,
            trial_inx = COMP{iD};
            switch g.PlotDashedCondition
                case 'min'
                    g.PlotDashedTime(iCELL,iD)=absmin(DashedTimes(trial_inx));
                case 'max'
                    g.PlotDashedTime(iCELL,iD)=absmax(DashedTimes(trial_inx));
                case 'mean'
                    g.PlotDashedTime(iCELL,iD)=nanmean(DashedTimes(trial_inx));
                case 'median'
                    g.PlotDashedTime(iCELL,iD)=nanmedian(DashedTimes(trial_inx));
            end
        end
    end
    
end % iC

% mean_npsth = squeeze(nanmean(shiftdim(NPSTH,1)));
% mean_npsth_se = squeeze(nanstd(shiftdim(NPSTH,1)))/sqrt(size(NPSTH,1)-1);

% dashtime = nanmean(unique(g.PlotDashedTime,'rows'));
% dashtime = nanmedian(unique(g.PlotDashedTime,'rows'));
dashtime = nanmax(unique(g.PlotDashedTime,'rows'));
% dashtime = nanmin(unique(g.PlotDashedTime,'rows'));
% dashtime = prctile(unique(g.PlotDashedTime,'rows'),100,1);

    for i=1:size(NPSTH,2)
        normNPSTH(i,:)    = nanmean(sq(NPSTH(:,i,:)));
        normNPSTHse(i,:) = nanstd(sq(NPSTH(:,i,:)))/sqrt(NumCells-1);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labelx=['Time-' g.TriggerName '(s)'];
tlimits = [g.window(1) g.window(2)];


figure(g.FigureNum)
if strcmp(g.Overlay,'off'),
    clf;
end
%subplot(211)
hold on;
for i=1:size(NPSTH,2)
    errorshade(time,normNPSTH(i,:),normNPSTHse(i,:),'Color',gcolor(i))
end
legend(Plegend)
pos2disp=restrict2(time,g.window(1),g.window(2));
MX = nanmax(nanmax(normNPSTH(:,pos2disp)))*1.15;
MN = nanmin(nanmin(normNPSTH(:,pos2disp)))*0.5;
alim=[g.window(1) g.window(2) MN MX];
axis(alim);
ylabel('Normalized rate');
xlabel(labelx);
setmyplot(gca);

figure(g.FigureNum+1)
if strcmp(g.Overlay,'off'),
    clf;
end
%subplot(211)
hold on;
lco = 1;
if ~isempty(g.PlotDashedEvent),
    for i = 1:size(NPSTH,2),
        if ~isnan(dashtime(i))
            posDASHED = nearest2(time,dashtime(i));
            if dashtime(i) < 0,
                p0(lco)=plot(time(posDASHED:end),normNPSTH(i,posDASHED:end),'Color',gcolor(i),'Tag',num2str(i));
                lco = lco+1;
                p0(lco)=plot(time(1:posDASHED),normNPSTH(i,1:posDASHED),g.DashedLineStyle,'Color',gcolor(i),'Tag',num2str(i));
                lco = lco+1;
            else
                p0(lco)=plot(time(1:posDASHED),normNPSTH(i,1:posDASHED),'Color',gcolor(i),'Tag',num2str(i));
                lco = lco+1;
                p0(lco)=plot(time(posDASHED:end),normNPSTH(i,posDASHED:end),g.DashedLineStyle,'Color',gcolor(i),'Tag',num2str(i));
                lco = lco+1;
            end
        else
            p0(lco)=plot(time,normNPSTH(i,:),'Color',gcolor(i),'Tag',num2str(i));
            lco = lco+1;
        end
    end
else
    for i=1:size(NPSTH,2)
        p1 = errorshade(time,normNPSTH(i,:),normNPSTHse(i,:),'Color',gcolor(i),'Tag',num2str(i));
        if  strcmp(g.PSTHstd,'off'),
            delete(p1(1));
            p0(lco) = p1(2);
            lco = lco+1;
        else
           disp('code this in')
        end
    end
end

legend(Plegend)
pos_disp=restrict2(time,g.window(1),g.window(2));
MX = max(max(normNPSTH(:,pos_disp)))*1.15;
MN = min(min(normNPSTH(:,pos_disp)))*0.5;
alim=[g.window(1) g.window(2) MN MX];
axis(alim);
ylabel('Normalized rate');
xlabel(labelx);
setmyplot(gca);
handles = p0;

% col={'g','r',colors6{3},colors6{6},'c','m','k',colors6{:},'b','y',colors6
