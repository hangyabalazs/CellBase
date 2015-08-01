function viewpopraster(cellids,partition,varargin)
%
%   VIEWPOP
%
%

if (nargin < 1)
	help viewpoppsth
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
TEfields = fieldnames(TE);

Npart = sum(ismember(lower(TEfields),lower(partition)));
if ismember('mixturediff',lower(partition))
    Npart=Npart+3;
end
if ismember('odorratio',lower(partition))
    Npart=Npart+6;
end

NPSTH=zeros(NumCells,Npart,length(time))*NaN;

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
end

alltrials = 1:size(ST.event_stimes{1},2);

stimes  = ST.event_stimes{trigger_pos}(alltrials);
windows = ST.event_windows{trigger_pos}(:,alltrials);

if strcmp(g.BurstPSTH,'on'),
    stimes = detect_bursts(stimes);
end

%%% MAKE THE MAIN RASTER
binraster = stimes2binraster(stimes,time,g.dt);

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
        NORMfactor =  mean(spsth(:));
    case 'median'
        NORMfactor =  median(spsth(:));
    case 'max'
        NORMfactor =  max(spsth(:));
    case 'min'
        NORMbaseline = nanmean(spsth(:));
        NORMfactor =  -min(spsth(:)-NORMbaseline);
    case 'zscore'
        NORMfactor =  std(spsth(:));
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

end % iC

% mean_npsth = squeeze(nanmean(shiftdim(NPSTH,1)));
% mean_npsth_se = squeeze(nanstd(shiftdim(NPSTH,1)))/sqrt(size(NPSTH,1)-1);
 
for i=1:size(NPSTH,2)
    normNPSTH(i,:)    = nanmean(sq(NPSTH(:,i,:)));
    normNPSTHse(i,:) = nanstd(sq(NPSTH(:,i,:)))/sqrt(NumCells-1);
end

% SORT
% this will only work correctly for now if there are no partitions


if ~isfield(g,'SortThreshold')
    g.SortThreshold =  prctile(NPSTH(:),66);
end

for iC=1:size(NPSTH,1)

    switch lower(g.Sort)
    case 'max'
        psth_max =  find(NPSTH(iC,1,:)==max(NPSTH(iC,1,:)),1);
        if ~isnan(psth_max)
            sort_ind(iC) = psth_max;
        else
            sort_ind(iC) = 1;
        end
    case 'onset'
      sort_ind(iC) = find(NPSTH(iC,1,:)>g.SortThreshold,1);
    otherwise
            sort_ind(iC) = iC;
    end %switch
end %iC

[junk reindex]=sort(sort_ind);
psth_raster = squeeze(NPSTH(:,1,:));
psth_raster = psth_raster(reindex,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labelx=['Time-' g.TriggerName '(s)'];
tlimits = [g.window(1) g.window(2)];



figure(g.FigureNum)
if strcmp(g.Overlay,'off'),   
clf;
end
subplot(311)
imagesc(time,1:size(psth_raster,1),psth_raster);
caxis([-max(abs(caxis)) max(abs(caxis))])
colorbar;
axis xy;
xlabel(labelx);
ylabel('Neuron #');
ax1=gca;
setmyplot(ax1);
subplot(313)
hold on;
for i=1:size(NPSTH,2)
    errorshade(time,normNPSTH(i,:),normNPSTHse(i,:),'Color',gcolor(i))
end
legend(Plegend)
pos_disp=restrict2(time,g.window(1),g.window(2));
MX = max(max(normNPSTH(:,pos_disp)))*1.15;
MN = min(min(normNPSTH(:,pos_disp)))*0.5;
alim=[g.window(1) g.window(2) MN MX];
%alim=[g.window(1) g.window(2) -0.7 2.0];
axis(alim);
ylabel('Normalized rate');
xlabel(labelx);
setmyplot(gca);

pos_ax1 = get(ax1,'Position');
set(ax1,'Position',[pos_ax1(1) pos_ax1(2)-pos_ax1(4)*1.3 pos_ax1(3) pos_ax1(4)*2.3])

% col={'g','r',colors6{3},colors6{6},'c','m','k',colors6{:},'b','y',colors6
