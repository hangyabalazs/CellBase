function viewraster(cellid,Partitions,varargin)
%
%   ViewRaster
%
%
if (nargin < 1)
	help viewraster
	return
end

% check if cellid is valid  
% if  validcellid(cellid,{'list'}) ~= 1
%     fprintf('%s is not valid.',cellid);
%     return
% end

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
        'ShowEvents',       {{'OdorPokeIn','OdorPokeOut','WaterPokeIn','WaterValveOn','WaterPokeOut'}};...
        'ShowEventsColors', {{'c','m','y','b','r'}};...
        'OdorPairID'        1;  ...
        'Normalization'     'max'; ...
        'NormalizationWindow'    []; ...
        'NormalizationTrials'    'all'; ...
        'window'            [-0.5 1]; ...
        'dt'                0.01; ...
        'sigma'             0.02; ...
        'plot'              'on'; ...
        'FigureNum'         1; ...
        'ValidTrials'       ''; ...
    };

[g, error] = parse_args(default_args,varargin{:});
g.Partitions = Partitions;

% test arguments for consistency
switch lower(g.plot)
    case { 'on', 'off' }, ;
    otherwise error('PLOT must be either on or off');
end;


%%-----------------------------------------------------------
%%  Preprocessing
%%-----------------------------------------------------------

margin = g.sigma*3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;  % time base array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cellid='N49_050803_4.1';

TE = loadcb(cellid,'Events');
SP = loadcb(cellid,'EVENTSPIKES');

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

%%% Could be put as an option
if g.OdorPairID == 0
    valid_trials = selecttrial(TE,sprintf('OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.ValidTrials));
else
    valid_trials = selecttrial(TE,sprintf('OdorPairID == %d & OdorConc == 100 & OdorPokeValid & WaterPokeValid %s',g.OdorPairID,g.ValidTrials));
end

[psth, spsth, spsth_se] = binraster2psth(binraster,g.dt,g.sigma,COMPTRIALS,valid_trials);
  
% partition2 = 'Outcome';
% g.partition = 'Stimulus';
% 
% [ind, part1, part2] = trialsort(TE,g.partition,partition2,g.TriggerEvent,g.SortEvent,valid_trials);
% trial_order = valid_trials(ind);
% partitions{1} = part1;
% partitions{2} = part2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% m=make3color([0 1 0],[0 0 1],[1 0 0]);
% junk=m(1:11:end,:);
% %colors6=mat2cell(junk,[1 1 1 1 1 1]);
% colors6 = num2cell(junk',[1 3]);
% 
% switch lower(g.partition)
%     case 'stimulus'
%                colors1 = colors6;
%     case 'direction'
%                colors1 = {'b','r'};
% end
% 
% colors2 = {'g','r'};
% 
% partition_colors{1} = colors1;
% partition_colors{2} = colors2;
% 
%ShowEvents = {'OdorPokeIn','OdorPokeOut','WaterPokeIn','WaterValveOn','WaterPokeOut'};

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
else
    eval(['sort_var=TE.' g.SortEvent ' - TE.' g.TriggerEvent ';']);
end


%plot_raster(fighandle,time,binraster,trial_order,EventTimes,tlimits,partitions,partition_colors,labelx,labely);

[mylabels, mycolors, mycolors2,mylinestyle] = makeColorsLabels(@defineLabelsColors_default,TAGS);
XLabel = ['Time - ' g.TriggerEvent];
YLabel = 'Rate (Hz)';
Num2Plot = NaN;

fhandle0 = plot_raster2(stimes,time,valid_trials,COMPTRIALS,mylabels,EventTimes,window_margin,ev_windows,sort_var,g,'Colors',{mycolors},'Colors2',{mycolors2},'NumTrials2Plot',Num2Plot);

plot_timecourse(time,spsth,spsth_se,g,'FigureNum',fhandle0,'Colors',{mycolors},'LineStyle',{mylinestyle},'Legend',{mylabels},'XLabel',XLabel,'YLabel',YLabel);


%---
% make fast? later

%---

% 
% figure(g.FigureNum+1)
% clf;
% col={[0 0 1],[1 0 0],colors6{3},colors6{6},[0 1 1],[1 1 0],[1 1 1],colors6{:},[0 1 0],[1 0 1],colors6{:},colors6{:}};
% subplot(221)
%    hold on;
%    for i=1:2
%      errorshade(time,spsth(i,:),spsth_se(i,:),'LineColor',col{i});
%    end
%    legend('Correct','Error',2);
%    pos_disp=restrict(time,g.window(1),g.window(2));
%    MX = (max(max(spsth(:,pos_disp))) +  max(max(spsth_se(:,pos_disp))))*1.1;
%    alim=[g.window(1) g.window(2) 0 MX];
%    axis(alim);
%    xlabel(labelx); ylabel('Firing rate');
% subplot(222)
%    hold on;
%    for i=3:4
%      errorshade(time,spsth(i,:),spsth_se(i,:),'LineColor',col{i});
%    end   
%    legend('Left','Right',2);
%    axis(alim);
%    xlabel(labelx); ylabel('Firing rate');
% subplot(223)
%    hold on;
%    for i=5:8
%      errorshade(time,spsth(i,:),spsth_se(i,:),'LineColor',col{i});
%    end   
%    legend('LC','LE','RC','RE',2);
%    axis(alim);
%    xlabel(labelx); ylabel('Firing rate');
% subplot(224)
%    hold on;
%    for i=9:size(spsth,1)
%      %errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
%      plot(time,spsth(i,:),'color',col{i});
%    end     
%    legend(num2str([1:6]'),2)
%    axis(alim)
%    xlabel(labelx); ylabel('Firing rate');
%    
% figure(g.FigureNum+2)
% clf;
% %col={'g','r',colors6{3},colors6{6},'c','c--','m','m--',colors6{:}};
% subplot(211)
%    hold on;
%    for i=1:2
%      errorshade(time,spsth(i,:),spsth_se(i,:),'LineColor',col{i});
%    end
%    legend('Correct','Error',2);
%    pos_disp=restrict(time,g.window(1),g.window(2));
%    MX = (max(max(spsth(:,pos_disp))) +  max(max(spsth_se(:,pos_disp))))*1.1;
%    alim=[g.window(1) g.window(2) 0 MX];
%    axis(alim);
%    xlabel(labelx); ylabel('Firing rate');
% subplot(212)
%    hold on;
%    for i=9:size(spsth,1)
%      %errorshade(time,spsth(i,:),spsth_se(i,:),col{i});
%      plot(time,spsth(i,:),'color',col{i},'LineWidth',2);
%    end     
%    legend(num2str([1:6]'),2)
%    axis(alim);
%    xlabel(labelx); ylabel('Firing rate');