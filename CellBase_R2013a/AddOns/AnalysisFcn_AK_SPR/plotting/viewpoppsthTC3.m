function viewpoppsthTC3(cellids,Partitions,varargin)
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
        'TriggerEvent'     'WaterPokeIn'; ...
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
        'TCstart'           -0.2;...
        'TCend'             0.3;...
        'NormalizationTC'   'max';...
        'NormCAF'           'linear';...
        'Bins'              8;...
        'MinPoints'         0;...
        'MinRate'           -1;...
    };


[par, error] = parse_args(default_args,varargin{:});
par.Partitions = Partitions;
par.Partitions_SAR ={'#OdorRatio:{[0] [32] [44] [56] [68] [100]}'};
par.RatePrctiles = [50]; %[33 66];
   
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

PartNum=zeros(1,length(par.Partitions));
for iUC=1:length(uniqsescells)
    TE = loadcb(uniqsescells{iUC},'Events');
    [TRIALS, ATAGS{iUC},R,N,PartNumI] = partition_trials(TE,par.Partitions);
    NumParts(iUC) = length(TRIALS);
    PartNum = max([PartNum; PartNumI]);
end

PTAGS = unique([ATAGS{:}]);
NumPartitions = length(PTAGS);
%-------------------------------------
for iUC=1:length(uniqsescells)
    TE = loadcb(uniqsescells{iUC},'Events');
    [TRIALS, ATAGS_SAR{iUC},R,N,PartNumI] = partition_trials(TE,par.Partitions_SAR);
end

PTAGS_SAR = unique([ATAGS_SAR{:}]);
NumPartitions_SAR = length(PTAGS_SAR);
%-------------------------------------------
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

ind_time = restrict(time,par.TCstart-eps,par.TCend+eps);
%xunc0=0:0.05:1;
%xunc0=0.066:0.066:1;
xunc0=0.0:0.07:1;
%
NumRateRanges = length(par.RatePrctiles)+1;
RATE_STIMBYRATE = nan(NumCells,NumPartitions_SAR,NumRateRanges);
ACCU_STIMBYRATE = nan(NumCells,NumPartitions_SAR,NumRateRanges);
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
%      case 'maxrate'
%          NORMfactor = max(EpochRate(union(trialsC,trialsE)));
     otherwise
        NORMfactor = 1;
end

%here is the key to figure out which parts were calculated...
posPARTS = match_list(TAGS,PTAGS);
 
NPSTH(iCELL,posPARTS,1:length(time)) = spsth/ NORMfactor;

if ~isempty(par.PlotDashedCondition)
    PlotDashedTimeV(iCELL) = eval(par.PlotDashedCondition);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Stimulus-Accuracy-By Rate (SAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[COMPTRIALS_SAR, TAGS_SAR] = partition_trials(TE,par.Partitions_SAR,valid_trials);

for iCT = 1:length(COMPTRIALS_SAR)
    trials_SAR = COMPTRIALS_SAR{iCT};
    if length(trials_SAR) > 20
        rate_sar_all = sum([binraster(trials_SAR,ind_time)]');        
        correct_sar_all = TE.Correct(trials_SAR);
        posNONZERO = find(rate_sar_all>0);
        rate_sar_all = rate_sar_all(posNONZERO);
        correct_sar_all = correct_sar_all(posNONZERO);
        if ~isempty(rate_sar_all)
            edgesPERC =prctile(rate_sar_all,par.RatePrctiles);
            if edgesPERC == 0
                edgesPERC = mean(rate_sar_all);
            end
            edges = [-Inf edgesPERC Inf];
            
            [n_SAR,bins_SAR]=histc(rate_sar_all,edges);
            if length(rate_sar_all) ~= sum(n_SAR)
                disp('PROBLEM: length(rate_sar_all) ~= sum(n)')
            end
            %        n_SAR
            for iRANGE = unique(bins_SAR)
                trials_in_range = find(bins_SAR == iRANGE);
                posNONZERO = find(rate_sar_all(trials_in_range)>-1);
                
                RATE_SAR(iCT,iRANGE)    = nanmean(rate_sar_all(trials_in_range(posNONZERO)));
                %%%            RATE_SAR(iCT,iRANGE)    = length(posNONZERO)/length(trials_in_range);
                %  RATE_SAR_SE(iCT,iRANGE) = nanstd(rate_sar_all(trials_in_range))/sqrt(length(trials_in_range));
                ACCU_SAR(iCT,iRANGE)    = nanmean(correct_sar_all(trials_in_range(posNONZERO)));        
                %  ACCU_SAR_SE(iCT,iRANGE) = nanstd(correct_sar_all(trials_in_range))/sqrt(length(trials_in_range));
                
            end %iRANGE
        else%isempty rate_sar_all
            RATE_SAR(iCT,1:2) = NaN;
            ACCU_SAR(iCT,1:2) = NaN;
        end %isempty rate_sar_all
    else %isempty TRIALS_SAR
         RATE_SAR(iCT,1:2) = NaN;
         ACCU_SAR(iCT,1:2) = NaN;
    end
end %iCT

posPARTS = match_list(TAGS_SAR,PTAGS_SAR);

NORMfactor_SAR = max(RATE_SAR(:));
%RATE_STIMBYRATE = nan(NumCells,NumPartitions_SAR,NumRateRanges);
%ACCU_STIMBYRATE = nan(NumCells,NumPartitions_SAR,NumRateRanges);

RATE_STIMBYRATE(iCELL, posPARTS, 1:NumRateRanges) = RATE_SAR./NORMfactor_SAR;
ACCU_STIMBYRATE(iCELL, posPARTS, 1:NumRateRanges) = ACCU_SAR;


%_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

EpochRate = nanmean([binraster(:,ind_time)]');
posCORR = intersect(find(TE.Correct == 1),valid_trials);
posERR  = intersect(find(TE.Error == 1), valid_trials);
%
[caf , xbins , nt] = histratio(EpochRate(posCORR),EpochRate(posERR),par.Bins,'fixed_bins','minpoints',par.MinPoints); 
%[caf , xbins] = histratio2(EpochRate(posCORR),EpochRate(posERR),par.Bins);
NORMcaf(iCELL,:)     = normalize(xbins,caf,xunc0,par.NormCAF);    

if nanmean(EpochRate(posERR)) > 0.5
    %'h'
end

end % iC


meanNORMcaf     = nanmean(NORMcaf);
medNORMcaf      = nanmedian(NORMcaf);
seNORMcaf       = nanstd(NORMcaf)/sqrt(NumCells-1);

ActNumCells = sum(~isnan(NORMcaf));

seNORMcaf2       = nanstd(NORMcaf)./sqrt(ActNumCells-1);

for i=1:size(NPSTH,2)
    normNPSTH(i,:)   = nanmean(sq(NPSTH(:,i,:)));
    normNPSTHse(i,:) = nanstd(sq(NPSTH(:,i,:)))/sqrt(NumCells-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%% Tuning Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_time = restrict(time,par.TCstart-eps,par.TCend+eps);

x=NPSTH(:,:,ind_time);
RATE =sq(nanmean(shiftdim(x,2)));
RATE2 =sq(max(shiftdim(x,2)));
N=size(RATE,2);
switch lower(par.NormalizationTC)
    case 'max'
            RATE=RATE./repmat(max(RATE')',1,N);
    case 'mean'
            RATE=RATE./repmat(nanmean(RATE')',1,N);        
    case 'meanpsthmax'
            RATE=RATE./repmat(nanmean(RATE2')',1,N);    
    case 'maxpsthmax'
            RATE=RATE./repmat(max(RATE2')',1,N);   
end


RATEmn = nanmean(RATE);
RATEse = nanstd(RATE)/sqrt(NumCells-1);
%--------------
% RATE_STIMBYRATE(iCELL, posPARTS, 1:NumRateRanges) = RATE_SAR./NORMfactor_SAR;
% ACCU_STIMBYRATE(iCELL, posPARTS, 1:NumRateRanges) = ACCU_SAR;
ACCU_SARmean = sq(nanmean(ACCU_STIMBYRATE));
ACCU_SARse   = [[nanstd(ACCU_STIMBYRATE(:,:,1));  nanstd(ACCU_STIMBYRATE(:,:,2))]./sqrt([sum(~isnan(ACCU_STIMBYRATE(:,:,1))); sum(~isnan(ACCU_STIMBYRATE(:,:,2)))])]';
RATE_SARmean = sq(nanmean(RATE_STIMBYRATE));
RATE_SARse   = [[nanstd(RATE_STIMBYRATE(:,:,1));  nanstd(RATE_STIMBYRATE(:,:,2))]./sqrt([sum(~isnan(RATE_STIMBYRATE(:,:,1))); sum(~isnan(RATE_STIMBYRATE(:,:,2)))])]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------

[xlabels, xvalues] = makeXLabels(@defineXLabels_default,PTAGS);


for i=1:length(xvalues)
     xcorrect(i) = length(findstr('Correct',PTAGS{i}));
 end
posCORR = find(xcorrect==1);
posERR  = find(xcorrect==0);

UXV=unique(xvalues);
xvalues2 = xvalues;
for i=1:length(UXV)
    xvalues2 = replace(xvalues2,UXV(i),i);
end   
    
[xv1C, ind1]=sort(xvalues(posCORR));
xv2C = xvalues2(posCORR(ind1));
posC = posCORR(ind1);

[xv1E, ind2]=sort(xvalues(posERR));
xv2E = xvalues2(posERR(ind2));
posE = posERR(ind2);
%--------------
% RATE_SAR
%-------------
[xlabels_SAR, xvalues_SAR] = makeXLabels(@defineXLabels_default,PTAGS_SAR);
[xv_SAR ind_SAR] = sort(xvalues_SAR);
 
%-----------------------
titlestr2 = sprintf('%d cells - Norm:%s; TCnorm:%s; Win:[%0.1f,%0.1f]',NumCells,par.Normalization,par.NormalizationTC,par.TCstart,par.TCend);

fhandle1=figure(par.FigureNum);
clf; 
subplot(221)
hold on;
LW=[1.5 2];
MSize=8;
YLabel='Norm rate';
XLabel='Mixture ratio (%)';
% 'XTickLabels',{mylabels},
% plot_tuning(posERR1,[RATEmn(index21); RATEse(index21)]',par,'FigureNum',fhandle1,'MarkerSize',MSize,'Line','off','LineWidth',LW,'Color','k',...
%             'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',xv2);
plot_tuning(xv2C,[RATEmn(posC); RATEse(posC)]',par,'FigureNum',fhandle1,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color','g',...
            'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',xv1C);        
plot_tuning(xv2E,[RATEmn(posE); RATEse(posE)]',par,'FigureNum',fhandle1,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color','r',...
            'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',xv1E);
axis([0.5 6.5 min(RATEmn)*0.8 max(RATEmn)*1.1])
suptitle(titlestr2);
axis([0.5 6.5 min(RATEmn)*0.8 max(RATEmn)*1.1])
setmyplot(gca);
subplot(222)
hold on;
e=errorshade(xunc0,meanNORMcaf,seNORMcaf2,'LineColor',[0.4 0.4 0.8]);
axis([0 1 0.5 0.95])
xlabel('Normalized rate');
ylabel('Accuracy');
setmyplot(gca);
%--------------------------------------------------------------------------
subplot(223)
hold on;
LW=[1.5 2];
MSize=4;
YLabel='Norm rate';
XLabel='Mixture ratio (%)';

plot_tuning(1:length(xv_SAR),[RATE_SARmean(ind_SAR,1) RATE_SARse(ind_SAR,1)],par,'FigureNum',fhandle1,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color','b',...
            'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',xv_SAR);  
errorbar(1:length(xv_SAR),RATE_SARmean(ind_SAR,2),RATE_SARse(ind_SAR,2),'t0-r')
axis([0.5 6.5 0 max(RATE_SARmean(:))*1.2])
setmyplot(gca);
subplot(224)
hold on;
LW=[1.5 2];
MSize=4;
YLabel='Accuracy';
XLabel='Mixture ratio (%)';
plot_tuning(1:length(xv_SAR),[ACCU_SARmean(ind_SAR,1) ACCU_SARse(ind_SAR,1)],par,'FigureNum',fhandle1,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color','b',...
            'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',xv_SAR);  
errorbar(1:length(xv_SAR),ACCU_SARmean(ind_SAR,2),ACCU_SARse(ind_SAR,2),'t0-r')
axis([0.5 6.5 0.5 1])
setmyplot(gca);

%------------------------------------------------------------------
%------------------------------------------------------------------
function ncaf = normalize(bins,caf,Xunc,method)

%normalize
if size(bins,2) > 1
    bins = bins/max(bins);
    pNaN = find(~isnan(caf));
    if length(pNaN) > 1 
        if strcmp(method,'linear')          
            %         pNaN = setdiff(pNaN,[1 length(caf)]);
            %         caf(pNaN) = mean([caf(pNaN-1); caf(pNaN+1)]);
            ncaf = interp1(bins(pNaN),caf(pNaN),Xunc,'linear');
        else
            ncaf = interp1(bins(pNaN),caf(pNaN),Xunc,'nearest');
        end
    else
         ncaf = Xunc*NaN;
    end
else
    ncaf = Xunc*NaN;
end