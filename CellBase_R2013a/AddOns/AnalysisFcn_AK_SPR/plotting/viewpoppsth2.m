function viewpoppsth2(cellids,Partitions,varargin)
%
%   VIEWPOPPSTH
%
%

if (nargin < 2)
	help viewpoppsth
	return
end
% check if cellid is valid  validcellid(cellid,{'list'}) ~= 1




% assign defaults if no arguments passed
default_args = { ...
        'Partitions'        'All'; ...
        'Compute',          'psth'; ... %roc
        'Transform',        'swap';...
        'Cumulative',       'no';...
        'Windowsize',       0.12;...
        'Windowshift',      0.04;...
        'TriggerEvent'     'HomeZoneOut'; ...
        'LastEvents',        '';...
        'OdorPairID'        1;  ...
        'Normalization'     'max'; ...
        'NormalizationWindow'    []; ...
        'NormalizationTrials'    'all'; ...
        'window'            [-0.5 1]; ...
        'dt'                0.01; ...
        'sigma'             0.02; ...
        'plot'              'on'; ...
        'FigureNum'         1; ...
        'ClearFig'          'on'; ...
        'ValidTrials'       ''; ...
        'PlotDashedCondition' '';...
        'PlotDashedTime'    NaN;...
    };

[par, error] = parse_args(default_args,varargin{:});
par.Partitions = Partitions;

% test arguments for consistency
switch lower(par.plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;


%%-----------------------------------------------------------
%%  Preprocessing
%%-----------------------------------------------------------

NumCells = length(cellids);
margin = par.sigma*3;     % add an extra margin to the windows
time = par.window(1)-margin:par.dt:par.window(2)+margin;  % time base array


%figure out how many partitions we have this is 
%complicated by the fact that not all sessions will 
%have the same conditions and hence the same number 
%of partitions.
uniqsescells = unique_session_cells(cellids);

for iUC=1:length(uniqsescells)
    TE = loadcb(uniqsescells{iUC},'Events');
    [TRIALS, ATAGS{iUC}] = partition_trials(TE,par.Partitions);
    NumParts(iUC) = length(TRIALS);
end

% x=[ATAGS{:}]
% unique(char(x{:}),'rows')

PTAGS = unique([ATAGS{:}]);
NumPartitions = length(PTAGS);

if strcmpi(par.Compute,'roc')
    if NumPartitions == 2
        %good
        COMPUTE_ROC = 1;
        PTAGS = {['D:' PTAGS{1} '-' PTAGS{2}]};
        NPSTH=nan(NumCells,1,length(time));  %just 1 value
    else
        error('2 partitions are required for ROC');
    end
else
    COMPUTE_ROC = 0;
    NPSTH=nan(NumCells,NumPartitions,length(time));
end

if ~isempty(par.NormalizationWindow)
    nwindow_ind = [nearest(time,par.NormalizationWindow(1)):nearest(time,par.NormalizationWindow(2))];
else
    nwindow_ind = [1:length(time)];
end

%%-----------------------------------------------------------
%%  Loop across cells
%%-----------------------------------------------------------
fprintf('%d cells to process.\n',NumCells)
for iCELL = 1:NumCells
 
print_progress(iCELL,round(NumCells/100),5);

cellid = cellids{iCELL};
TE = loadcb(cellid,'Events');
ST = loadcb(cellid,'EVENTSPIKES');

trigger_pos = findcellstr(ST.events(:,1),par.TriggerEvent);

if (trigger_pos == 0)
  error('Trigger variable not found');
end


alltrials = 1:size(ST.event_stimes{1},2);
stimes  = ST.event_stimes{trigger_pos}(alltrials);;
windows = ST.event_windows{trigger_pos}(:,alltrials);



if ~iscellstr(par.LastEvents) & (strcmpi(par.LastEvents,'none') | isempty(par.LastEvents)) 
    window_margin = ST.events{trigger_pos,4};
    ev_windows = ST.event_windows{trigger_pos};          
else
    window_margin = [par.window(1)-2*par.dt 0];
    ev_windows = get_last_evtime(TE,par.TriggerEvent,par.LastEvents);
    % not very elegant, but we need the time before TriggerEvent
    ev_windows  =  ev_windows + repmat([time(1) 0],size(ev_windows,1),1); 
end

          

%%% MAKE THE MAIN RASTER
binraster = stimes2binraster(stimes,time,par.dt,ev_windows,window_margin);

[COMPTRIALS, TAGS] = partition_trials(TE,par.Partitions);

% if unique(TE.OdorRatio) == 4
%     %medium -> diff
%      junk = TAGS{3};
%      TAGS{3} = TAGS{4};
%      TAGS{4} = junk;
% end
%  

%%% Could be put as an option
if par.OdorPairID == 0
    valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',par.ValidTrials));
else
    valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',par.OdorPairID,par.ValidTrials));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(valid_trials) > 1
    if COMPUTE_ROC
        spsth = binraster2selectivitytime(binraster,COMPTRIALS,valid_trials,'Cumulative', ... 
            par.Cumulative,'Windowsize',par.Windowsize,'Windowshift',par.Windowshift,'Transform',par.Transform,'FigureNum',2);
        TAGS = sort(TAGS);
        TAGS = {['D:' TAGS{1} '-' TAGS{2}]};  %new tag required
    else
        [psth, spsth, spsth_se] = binraster2psth(binraster,par.dt,par.sigma,COMPTRIALS,valid_trials);
        
    end
else
    spsth    = nan(length(TAGS),length(time));
    spsth_se = spsth;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmpi(par.NormalizationTrials,'all')
    ALLTRIALS{1} = 1:length(TE.TrialStart);
    [psth, norm_spsth, spsth_se] = binraster2psth(binraster,par.dt,par.sigma,ALLTRIALS,valid_trials);
    norm_spsth = norm_spsth(nwindow_ind);
else
    norm_spsth = spsth(:,nwindow_ind);
end

switch lower(par.Normalization)
    case 'none'
        NORMfactor = 1;
    case 'mean'
        NORMfactor =  nanmean(norm_spsth(:));
    case {'median', 'medi'}
        NORMfactor =  nanmedian(norm_spsth(:));
    case 'max'
        NORMfactor =  max(norm_spsth(:));
    case 'perc90'   
         NORMfactor = prctile(norm_spsth(:),90);
%     case 'maxrate'
%         NORMfactor = max(EpochRate(union(trialsC,trialsE)));
    otherwise
        NORMfactor = 1;
end

%here is the key to figure out which parts were calculated...
posPARTS = match_list(TAGS,PTAGS);
 
NPSTH(iCELL,posPARTS,1:length(time)) = spsth/ NORMfactor;

if ~isempty(par.PlotDashedCondition)
        PlotDashedTimeV(iCELL) = eval(par.PlotDashedCondition);
end
% %  S{iCELL} =  unique(TE.Stimulus);
% %  OR{iCELL} =  unique(TE.OdorRatio); 
end % iC

% mean_npsth = squeeze(nanmean(shiftdim(NPSTH,1)));
% mean_npsth_se = squeeze(nanstd(shiftdim(NPSTH,1)))/sqrt(size(NPSTH,1)-1);
 
for i=1:size(NPSTH,2)
    normNPSTH(i,:)   = nanmean(sq(NPSTH(:,i,:)));
    normNPSTHse(i,:) = nanstd(sq(NPSTH(:,i,:)))/sqrt(NumCells-1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if iscell(par.LastEvents)
    LastEvents = char([par.LastEvents{:}]);
else
    LastEvents = '';
end

if isempty(par.NormalizationWindow)
    par.NormalizationWindow =par.window;
end

if ~isempty(par.PlotDashedCondition)
    par.PlotDashedTime = min(PlotDashedTimeV);
end
   


titlestr = sprintf('%d cells - Norm:%s=[%1.2f %1.2f] %s; LastEvs:%s',NumCells,par.Normalization,par.NormalizationWindow(1),par.NormalizationWindow(2),...
                    par.NormalizationTrials, LastEvents);
        
[mylabels, mycolors] = makeColorsLabels(@defineLabelsColors_default,PTAGS);
[mylabels,ind] = sort(mylabels);      
mycolors = mycolors(ind);
labelx=['Time - ' par.TriggerEvent ' (s)'];

switch lower(par.Normalization)
    case 'none'
        labely = 'Firing rate'; 
    otherwise
        labely = 'Normalized response';
end


f0 = plot_timecourse(time,normNPSTH(ind,:), normNPSTHse(ind,:),par,'XLabel',labelx,'YLabel',labely,'Legend',{mylabels},'Colors',{mycolors},'PSTHstd','on','TitleStr',titlestr);

% SEL=sq(NPSTH);
% Thres=[-2 2];
% par.NormalizationWindow = [-0.2 0.5];
% par.PlotAverage='on';
% f0=plot_selectivity(time,SEL,par,'Thres',Thres)
