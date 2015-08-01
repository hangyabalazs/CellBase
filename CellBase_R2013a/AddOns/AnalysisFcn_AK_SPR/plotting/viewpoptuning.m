function [NVALUES, CP] = viewpoptuning(cellids,Partitions,varargin)
%
%   VIEWPOPTUNING
%
%



% assign defaults if no arguments passed
default_args = { ...
        'EpochRate'        'AfterChoice'; ...
        'Analysis'          'nanmean(TE.Correct(trials))';...
%         'TriggerEvent'     'WaterPokeIn'; ...
%         'LastEvents',        '';...
        'OdorPairID'        1;  ...
        'Normalization'     'max'; ...
%         'NormalizationWindow'    []; ...
        'NormalizationTrials'    ''; ...   %'all'
%         'window'            [-0.5 1]; ...
        'Plot'              'on'; ...
        'FigureNum'         1; ...
        'ValidTrials'       ''; ...
    };

[g, error] = parse_args(default_args,varargin{:});
g.Partitions = Partitions;

% test arguments for consistency
switch lower(g.Plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumCells = length(cellids);

%figure out how many partitions we have this is 
%complicated by the fact that not all sessions will 
%have the same conditions and hence the same number 
%of partitions.
uniqsescells = unique_session_cells(cellids);

for iUC=1:length(uniqsescells)
    TE = loadcb(uniqsescells{iUC},'Events');
   % [TRIALS, ATAGS{iUC}, RTAGS, NUM_TAGS{iUC}] = partition_trials(TE,g.Partitions);
    [TRIALS, ATAGS{iUC}] = partition_trials(TE,g.Partitions);
    NumParts(iUC) = length(TRIALS);
end

PTAGS = unique([ATAGS{:}]);
%%NTAGS = unique([NUM_TAGS{:}]); -- Not sure what this is for
NumPartitions = length(PTAGS);

if ~isempty(g.Analysis)
    NumVals = 2;
else
    NumVals = 1;
end

NVALUES = nan(NumCells,NumPartitions,NumVals);

%NVALUES = NumCell x NumPartitions x [EpochRate Analysis]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%d cells to process.\n',NumCells)
for iCELL = 1:NumCells
 
print_progress(iCELL,round(NumCells/100),5);

cellid = cellids{iCELL};
TE = loadcb(cellid,'Events');
ST = loadcb(cellid,'EVENTSPIKES');

epoch_pos = findcellstr(ST.epochs(:,1),g.EpochRate);
if (epoch_pos == 0)
  error('Trigger variable not found');
end

EpochRate = ST.epoch_rates{epoch_pos};

%%%-----trialsC = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & Correct ',OdorPairID));
  
alltrials = 1:size(ST.event_stimes{1},2);
%%stimes  = ST.event_stimes{trigger_pos}(alltrials);;
%%windows = ST.event_windows{trigger_pos}(:,alltrials);

%[COMPTRIALS, TAGS, RTAGS, NUM_TAGS] = partition_trials(TE,g.Partitions);
[COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);

%%% Could be put as an option
if g.OdorPairID == 0
    valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.ValidTrials));
else
    valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.OdorPairID,g.ValidTrials));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xPROP = 1 x NumPartitions x [EpochRate Analysis]

[values, values_se, accu]= compare_trials(EpochRate,TE,g.Analysis,COMPTRIALS,valid_trials);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pn=find(~isnan(accu));
warning off MATLAB:divideByZero
[cc,p]=corrcoef(accu(pn),EpochRate(pn));
warning on MATLAB:divideByZero
CP(iCELL,:) = [cc(2,1) p(2,1)];
%% Hmmm...
% if length(NUM_TAGS) == 6 & CP(iCELL,2) < 0.02
%    % error('x');
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(g.NormalizationTrials,'all')
     ALLTRIALS{1} = 1:length(TE.TrialStart);
     [norm_vals, values_se]= compare_trials(EpochRate,TE,g.Analysis,ALLTRIALS,valid_trials);
else
     norm_vals = values;
end

switch lower(g.Normalization)
    case 'none'
        NORMfactor = 1;
    case 'mean'
        NORMfactor =  nanmean(norm_vals(:));
    case {'median','medi'}
        NORMfactor =  nanmedian(norm_vals(:));
    case 'max'
        NORMfactor =  max(norm_vals(:));
    case 'perc90'   
         NORMfactor = prctile(norm_vals(:),90);
    case 'maxrate'
         NORMfactor = max(EpochRate);
    case 'meanrate'
         NORMfactor = nanmean(EpochRate);
    case 'perc90rate'
         NORMfactor = max(prctile(EpochRate,98),nanmean(EpochRate));
    otherwise
         NORMfactor = 1;
end

%here is the key to figure out which parts were calculated...
posPARTS = match_list(TAGS,PTAGS);

% value(:,2) is accuracy, no need to normalize
NVALUES(iCELL,posPARTS,1:NumVals) = [values(:,1) / NORMfactor  values(:,2)];


end % iC

for i=1:size(NVALUES,2)
    normNVALUES(i,:)   = nanmean(sq(NVALUES(:,i,:)));
    normNVALUESse(i,:) = nanstd(sq(NVALUES(:,i,:)))/sqrt(NumCells-1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(g.Plot,'on')
    
[mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_default,PTAGS);

[mylabels,ind] = sort(mylabels);
for i=1:length(PTAGS)
    [s,rem]=strtok(PTAGS{i},'=');
    [rem2,rem3]=strtok(rem,'|&');
    if ~isempty(rem2)
        n_val(i) = str2num(rem2(2:end));
    else
        n_val(i) = NaN;
    end
end

if ~sum(isnan(n_val))
    [x_val ind]=sort(n_val);
    nope = find(x_val ~= 50);
    x_val = x_val(nope);
    ind   = ind(nope);
else
    x_val = NaN;
end

mNVALUES = NVALUES(:,ind,:);
titlestr = sprintf('%d cells; Norm:%s - %s',NumCells,g.Normalization,g.EpochRate);
        
        


fhandle1=figure(g.FigureNum+1);
clf;
shandle1=subplot(221);
shandle2=subplot(222);
shandle3=subplot(223);
shandle4=subplot(224);
%
LW=[1.5 2];
MSize=8;
YLabel='Norm rate';
XLabel='Mixture ratio (%)';
% 'XTickLabels',{mylabels},
mNVALUES = [normNVALUES(ind,1) normNVALUESse(ind,1)];

plot_tuning(1:length(x_val),mNVALUES,g,'FigureNum',shandle1,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color',{mycolors},...
            'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',x_val);
%axis square
%
%----------------
subplot(223)
e0 =errorbar(1:length(x_val),normNVALUES(ind,2),normNVALUESse(ind,2),'t0');
yl=xlabel('Accuracy (%)');
xl=xlabel('Mixture ratio (%)');
XMIN=min(min([normNVALUES(ind,2)-normNVALUESse(ind,2)]))*0.95;
axis([0.5 length(x_val)+0.5 XMIN 1]);
set(gca,'XTickLabel',x_val);
set([e0],'Color','k');
set(e0(1),'LineWidth',1);
set(e0(2),'LineWidth',2);
box off;
setmyplot(gca);
%----------------------------
XLabel='Accuracy';
YLabel='';
plot_tuning(normNVALUES(ind,2),mNVALUES,g,'FigureNum',shandle2,'MarkerSize',MSize,'Line','off','LineWidth',LW,'Limits',[],'Color',{mycolors},...
            'XLabel',XLabel,'YLabel',YLabel);
drawnow
suptitle(titlestr);

% titstr=sprintf('StimAlign %s; Center %s; Normalize %s',ONOFF{g.StimAlign+1},ONOFF{g.StimCenter+1},g.StimNormalization);
% title(titstr);
% setmyplot(gca);
end