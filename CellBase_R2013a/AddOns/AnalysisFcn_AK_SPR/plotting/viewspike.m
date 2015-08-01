function viewspike(cellid,varargin)
% viewspike(cellid,varargin)
% plots for a given cell, the mean waveform (light triggered and
% non-triggered)
% Cluster quality
% Waveform properties

% inputs:
% cellid: (required)

if (nargin < 1)
    help viewspike
    return
end

default_args={...
    'window',               [-0.5 1];...
    'dt',                   0.01;...
    'sigma',                0.02;...
    'FigureNum',            1;...
    'TriggerEvent',         'PulseOn';...
    'SortEvent',            'PulseOff';...
    'filetype',             'stim';...
    'ShowEvents',           {{'PulseOff'}};...
    'ShowEventsColors',     {[1 0 0]};...
    'Num2Plot',             100;...
    ...
    'PlotDashedCondition',  '';...
    'LastEvents',           '';...
    'Partitions',           'all';...
    'PrintCellID',          'on';...
    };
[g,error]=parse_args(default_args,varargin{:});


% check if cellid is valid
if  validcellid(cellid,{'list'}) ~= 1
    fprintf('%s is not valid.',cellid);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cellid='N49_050803_4.1';
% StimEvent, TrialEvent mismatch
if ~isfield(TE,g.TriggerEvent),
    error('TriggerEvent mismatch: supply correct Events structure')
end

trigger_pos = findcellstr(SP.events(:,1),g.TriggerEvent);

if (trigger_pos == 0)
    error('Trigger variable not found');
end

alltrials = 1:size(SP.event_stimes{1},2);
stimes  = SP.event_stimes{trigger_pos}(alltrials);
windows = SP.event_windows{trigger_pos}(:,alltrials);


if ~iscellstr(g.LastEvents) & (strcmpi(g.LastEvents,'none') | isempty(g.LastEvents))
    window_margin = SP.events{trigger_pos,4};
    ev_windows = SP.event_windows{trigger_pos};
else
    window_margin = [g.window(1)-2*g.dt 0];
    ev_windows = get_last_evtime(TE,g.TriggerEvent,g.LastEvents);
end

%%% MAKE THE MAIN RASTER
binraster = stimes2binraster(stimes,time,g.dt,ev_windows,window_margin);

NUM_TRIALS = length(alltrials);

margin = g.sigma*3; % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array

%%% MAKE THE MAIN RASTER
binraster = stimes2binraster(stimes,time,g.dt);

[COMPTRIALS, TAGS] = partition_trials(TE,g.Partitions);
valid_trials=1:length(~isnan(getfield(TE,g.TriggerEvent)));
%%% Could be put as an option
% if g.OdorPairID == 0
%     valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.ValidTrials));
% else
%     valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.OdorPairID,g.ValidTrials));
% end

[psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);

%--------------------------------------------------------

EventTimes = trialevents2relativetime(TE,g.TriggerEvent,g.ShowEvents);

labelx=['Time-' g.TriggerEvent];
labely='Trials';
fighandle = g.FigureNum;
tlimits = [g.window(1) g.window(2)];

if iscellstr(g.SortEvent)
    for iS=1:length(g.SortEvent)
        eval(['sort_var(iS,:)=TE.' g.SortEvent{iS} ' - TE.' g.TriggerEvent ';']);
    end
    sort_var = min(sort_var);
elseif ~isempty(g.SortEvent)
    eval(['sort_var=TE.' g.SortEvent ' - TE.' g.TriggerEvent ';']);
else
    sort_var = NaN;
end


%plot_raster(fighandle,time,binraster,trial_order,EventTimes,tlimits,partitions,partition_colors,labelx,labely);

[mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_default,TAGS);
XLabel = ['Time - ' g.TriggerEvent];
YLabel = 'Rate (Hz)';


%-----------------------------
%   Raster + psth
%-----------------------------
fhandle0 = plot_raster2(stimes,time,valid_trials,COMPTRIALS,mylabels,EventTimes,window_margin,ev_windows,sort_var,g,'Colors',{mycolors},'Colors2',{mycolors2},'NumTrials2Plot',g.Num2Plot);

if ~isempty(g.PlotDashedCondition)
    %junk = valid_trials;
    
    g.PlotDashedTime = eval(g.PlotDashedCondition);
end
plot_timecourse(time,spsth,spsth_se,g,'FigureNum',fhandle0,'Colors',{mycolors},'LineStyle',{mylinestyle},'Legend',{mylabels},'XLabel',XLabel,'YLabel',YLabel);

%-----------------------------
%   Tuning + caf
%-----------------------------
% Partitions = {'#OdorRatio'};
% [COMPTRIALS1, TAGS, RTAGS, NUM_TAGS] = partition_trials(TE,Partitions);
% [mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_default,TAGS);
% 
% Accuracy = get_tuning(TE.Correct,{'nanmean(data)';'length(data)'},COMPTRIALS1,valid_trials);
% 
% [a,v]=binostat(1,Accuracy(:,1));
% Accuracy(:,2) = sqrt(v)./sqrt(Accuracy(:,2)-1);
% %-------------------------
% 
% 
% trigger_pos = findcellstr(SP.epochs(:,1),g.EpochName);
% if (trigger_pos == 0)
%     error('Epoch not found');
% end
% epoch_rate = SP.epoch_rates{trigger_pos};
% MeanRate = get_tuning(epoch_rate,{'nanmean(data)';'nanstd(data)/sqrt(length(data))'},COMPTRIALS1,valid_trials);
% %-----------------------------
% %   Stimulus-Accuracy-by Rate (SAR)
% %-----------------------------
% 
% for iCT = 1:length(COMPTRIALS1)
%     rate_sar_all = epoch_rate(COMPTRIALS1{iCT});
%     correct_sar_all = TE.Correct(COMPTRIALS1{iCT});
%     %    prctile_groups = [33 66];
%     %    prctile_groups = [66];
%     edges = [-Inf prctile(rate_sar_all,g.RatePrctiles) Inf];
%     [n_SAR,bins_SAR]=histc(rate_sar_all,edges);
%     if length(rate_sar_all) ~= sum(n_SAR)
%         disp('PROBLEM: length(rate_sar_all) ~= sum(n)')
%     end
%     for iRANGE = unique(bins_SAR)
%         trials_in_range = find(bins_SAR == iRANGE);
%         posNONZERO = find(rate_sar_all(trials_in_range)>g.MinRate);
%         RATE_SAR(iCT,iRANGE)    = mean(rate_sar_all(trials_in_range(posNONZERO)));
%         RATE_SAR_SE(iCT,iRANGE) = nanstd(rate_sar_all(trials_in_range(posNONZERO)))/sqrt(length(trials_in_range(posNONZERO)));
%         ACCU_SAR(iCT,iRANGE)    = mean(correct_sar_all(trials_in_range(posNONZERO)));
%         ACCU_SAR_SE(iCT,iRANGE) = nanstd(correct_sar_all(trials_in_range(posNONZERO)))/sqrt(length(trials_in_range(posNONZERO)));
%     end %iRANGE
% end %iCT
% 
% %-------------------------------------------
% if strcmpi(g.PrintCellID,'on')
%     fstamp(cellid,80);
% end
% 
% fhandle1=figure(g.FigureNum+1);
% clf;
% shandle1=subplot(221);
% shandle2=subplot(222);
% shandle3=subplot(223);
% shandle4=subplot(224);
% %
% LW=[1.5 2];
% MSize=0;
% %
% YLabel='Firing rate';
% XLabel='Stimulus';
% % 'XTickLabels',{mylabels},
% plot_tuning(1:length(NUM_TAGS),[RATE_SAR(:,2) RATE_SAR_SE(:,2)],g,'FigureNum',shandle1,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color',{mycolors},...
%     'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',NUM_TAGS,'TitleStr',g.EpochName);
% errorbar(1:length(NUM_TAGS),RATE_SAR(:,1),RATE_SAR_SE(:,1),'t0-b')
% axis([0.5 6.5 0 round(max(RATE_SAR(:))*1.2)]);
% %axis square
% %
% 
% XLabel='STIMULUS';
% YLabel='Accuracy';
% plot_tuning(1:length(NUM_TAGS),[ACCU_SAR(:,2) ACCU_SAR_SE(:,2)],g,'FigureNum',shandle2,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Color',{mycolors},...
%     'XLabel',XLabel,'YLabel',YLabel,'XTickLabels',NUM_TAGS,'TitleStr',g.EpochName);
% errorbar(1:length(NUM_TAGS),ACCU_SAR(:,1),ACCU_SAR_SE(:,1),'t0-b')
% axis([0.5 6.5 0.5 1]);
% %axis square
% %-----------------------------------
% Partitions = {'#OdorRatio:{0 20 32 44 50 56 68 80 100} & Correct','#OdorRatio:{0 20 32 44 50 56 68 80 100} & Error'};
% [COMPTRIALS2, TAGS2] = partition_trials(TE,Partitions);
% STIMULI = [0 20 32 44 50 56 68 80 100];
% MeanRateCE = get_tuning(epoch_rate,{'nanmean(data)';'nanstd(data)/sqrt(length(data))'},COMPTRIALS2,valid_trials);
% 
% MXX=max([MeanRateCE(:,1)+MeanRateCE(:,2)]);
% alim = [NaN NaN 0 MXX];
% ind=[];
% for iC=1:length(COMPTRIALS2)
%     if ~isempty(COMPTRIALS2{iC})
%         ind=[ind iC];
%     end
% end
% lenSTIM = length(STIMULI);
% posC=intersect(1:lenSTIM,ind);
% posE=intersect(lenSTIM+1:lenSTIM*2,ind);
% 
% plot_tuning(1:length(posC),MeanRateCE(posC+lenSTIM,:),g,'FigureNum',shandle3,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Limits',alim,'Color',[0.8 0.2 0.1]);
% plot_tuning(1:length(posC),MeanRateCE(posC,:),g,'FigureNum',shandle3,'MarkerSize',MSize,'Line','on','LineWidth',LW,'Limits',alim,'Color',[0.2 0.8 0.1],...
%     'XTickLabels',STIMULI,'XLabel','Stimulus','YLabel','Firing rate','TitleStr','C/E');
% 
% 
% 
% %---------------
% 
% NumBins = 6;
% %NumBins = 11:6.5:50;
% NumBins = 0.5:4:35;
% MinPoints = 5;
% Method = 'fixed_bins';
% [cafA, cafA_se, xbinsA, nt] = get_caf(epoch_rate,TE,{'Correct','Error'},valid_trials,NumBins,MinPoints,Method);
% [cafL, cafL_se, xbinsL] = get_caf(epoch_rate,TE,{'LeftCorrect','LeftError'},valid_trials,NumBins,MinPoints,Method);
% [cafR, cafR_se, xbinsR] = get_caf(epoch_rate,TE,{'RightCorrect','RightError'},valid_trials,NumBins,MinPoints,Method);
% 
% if ~min(isnan(xbinsA)) & ~min(isnan(cafA))
%     focusfigure(shandle4);
%     %errorshade(xbinsA,cafA,cafA_se,'LineColor',[0.2 0.2 0.8],'LineWidth',2.5);
%     %pcaf=plot(xbinsA,cafA,'o');
%     %set(pcaf,'MarkerEdgeColor','none','MarkerFaceColor',[0.2 0.2 0.8],'MarkerSize',MSize);
%     
%     errorshade(xbinsL,cafL,cafL_se,'LineColor',[0.3 0.2 0.6],'LineWidth',2);
%     errorshade(xbinsR,cafR,cafR_se,'LineColor',[0.3 0.6 0.2],'LineWidth',2);
%     pL=plot(xbinsL,cafL,'o');
%     pR=plot(xbinsR,cafR,'o');
%     MSize = 2;
%     set(pR,'MarkerEdgeColor','none','MarkerFaceColor',[0.3 0.6 0.2],'MarkerSize',MSize);
%     set(pL,'MarkerEdgeColor','none','MarkerFaceColor',[0.3 0.2 0.6],'MarkerSize',MSize);
%     
%     plot(xbinsA([1 end]),[0.5 0.5],'k--');
%     %set(pL,'Color','c','LineWidth',2);
%     %set(pR,'Color','y','LineWidth',2);
%     l(1)=xlabel('Firing rate');
%     l(2)=ylabel('Accuracy');
%     l(3)=title('CAF (all/left/right)');
%     Y=[cafA-cafA_se cafA+cafA_se cafL cafR];
%     X=[xbinsA xbinsL xbinsR];
%     xMN=min(X(:)); xMX=max(X(:));
%     yMN=min(Y(:)); yMX=max(Y(:));
%     xDX=(xMX-xMN)*0.05; yDX=(yMX-yMN)*0.1;
%     yDX=max(2*eps,yDX); xDX=max(2*eps,xDX);
%     axis([xMN-xDX xMX+xDX yMN-yDX yMX+yDX]);
%     setmyplot(gca,l);
% end
% 

%---------------------
if strcmpi(g.PrintCellID,'on')
    fstamp(cellid,80);
end