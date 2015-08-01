function fhandle0 = plot_raster2a(binraster,stimes,time,valid_trials,COMPTRIALS,TAGS,EventTimes,window_margin,ev_windows,sort_var,params,varargin)
%PLOT_RASTER2A   Plots raster.
%   H = PLOT_RASTER2A(STIMES,TIME,VALID_TRIALS,COMPTRIALS,TAGS,EVENTTIMES,WINDOW_MARGIN,EV_WINDOWS,SORT_VAR,PARAM,VALUE)
%   plots a raster on an axis/subplot within a figure. Input arguments:
%       STIMES: spike times
%       VALID_TRIALS: trials ton include in the plot
%       COMPTRIALS: trial partitions, see PARTITION_TRIALS
%       TAGS: tags corresponding to trial partitions, see PARTITION_TRIALS
%       EVENTTIMES: times of events to show
%       WINDOW_MARGIN: a margin applied to all windows
%       EV_WINDOWS: time frames to use around the trigger stimuli
%       SORT_VAR: trials are sorted according to this variable
%       additional PARAM VALUE pairs can be passed to the function. See the
%       default parameter settings in the code to get a list of optional
%       parameters.
%   Handle of the output figure is returned as output argument.
%
%   See also VIEWCELL2B and PARTITION_TRIALS.

%   Edit log: BH 6/27/11

% Set figure rendering
set(gcf,'renderer','painters')

% Set default parameters
default_args = params;
default_args.TitleStr = '';
default_args.XLabel = 'Time';
default_args.PSTH = 'on';
default_args.ShowRasterTAGS=1;
default_args.NumTrials2Plot = 20;
default_args.RasterIDbarDX  = 0.01;
default_args.RasterIDbarWidth = 6;
default_args.RasterPanelDistance = 0.08;
default_args.PsthPanelDistance=0.04;
default_args.RasterPlotNumtrials = 1;
default_args.EventMarkerWidth = 0.008;
default_args.EventMarkerHeight = 2;
default_args.SpikeWidth = 0.5;
default_args.ForceSort = 1;    % sort even in the presence of NaN values (good for Previous and Next Events where either the first or last value is a NaN
default_args.PlotZeroLine = 'off';
default_args.PlotZeroLineColor = [0 0 0];
default_args.PrintCellID='off';
default_args.NumPSTHPlots=NaN;
[par,error] = parse_args(default_args,params);   % overwrite defaults with passed arg structure
[par,error] = parse_args(par,varargin{:});

% Event windows
NumTrials = length(stimes);
if strcmp(par.NumTrials2Plot,'all'),
    par.NumTrials2Plot = NumTrials;
end
if ~exist('window_margin','var')
   window_margin = [0 0];
end
if exist('ev_windows','var')
    if size(ev_windows,1) ~= NumTrials
        ev_windows = ev_windows';
    end
    if ~iscell(ev_windows)   % deals with lick-aligned rasters - BH
        ev_windows = ev_windows + repmat(window_margin,NumTrials,1);
    end
end
det = par.EventMarkerWidth / 2;
dnt = par.EventMarkerHeight / 2;

% Initialize figure
NumParts = length(COMPTRIALS);
if ~isnan(par.NumPSTHPlots),
    NumPsthPlots = par.NumPSTHPlots;
else
    NumPsthPlots = 1;
end

if strcmpi(par.PSTH,'on')
    NumPanels = NumParts + NumPsthPlots;
else
    NumPanels = NumParts;
end

if ishandle(par.FigureNum),
    if strcmp(get(par.FigureNum,'Type'),'figure')
        clf;
        fhandle0 = axes;
    elseif strcmp(get(par.FigureNum,'Type'),'axes')
        fhandle0 = par.FigureNum;
        axes(fhandle0)
        cla;
    end
    ax = splitaxes(fhandle0,NumPanels);
else
    fhandle0 = focusfigure(par.FigureNum,'create');
    fhandle0 = axes;
    ax = splitaxes(fhandle0,NumPanels);
end

% Plot
iAX = 1;
for iCOMP = 1:NumParts
    
    % Get trials to plot
    ind = intersect(COMPTRIALS{iCOMP},valid_trials);    % use only valid ones
    if par.NumTrials2Plot < NumTrials
        rndind = randperm(length(ind));    % permute trials only if not all of them are plotted
    else
        rndind = 1:length(ind);
    end
    ind = ind(rndind(1:min(length(ind),par.NumTrials2Plot)));
    NumTrials(iCOMP) = length(ind);
        
    % Reindex
    if par.ForceSort == 1
        [i,sind] = sort(sort_var(ind));   % sort by the sort variable
        ind = ind(sind);
    else    % do it only if there are no NaN values
        if ~isnan(sort_var)    % reindex
            [i,sind] = sort(sort_var(ind));    % sort by the sort variable
            ind = ind(sind);
        end
    end
    
    % Plot spikes
    if NumTrials(iCOMP) > 0
        if ~iscell(ev_windows)
            nt = 1:NumTrials(iCOMP);
        else
            stimes2 = stimes(ind);
            numallicks = sum(cellfun(@length,stimes2));
            nt = 1:numallicks;
        end
        tbinraster = binraster(ind,:);
        S = subplot(ax(iAX));
        rasterplot(tbinraster,time,S)
        
        % Plot ShowEvents
        seti = size(EventTimes,1);
        p = zeros(1,seti);
        for iE = 1:seti
            if ~iscell(ev_windows)
                if ~iscell(EventTimes)
                    et = EventTimes(iE,ind);
                    et(et<time(1)) = NaN;
                    et(et>time(end)) = NaN;
                    if length(find(et==0)) == NumTrials(iCOMP)
                        line([0 0],[1 NumTrials(iCOMP)],'Color',par.ShowEventsColors{iE},'LineWidth',1);
                    else
                        p(iE) = patch([et-det; et-det; et+det; et+det],[nt+dnt; nt-dnt; nt-dnt; nt+dnt],'k');
                        try
                            set(p(iE),'FaceColor',par.ShowEventsColors{iE},'EdgeColor','none');
%                             uistack(p(iE),'bottom')
                        catch
                            set(p(iE),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.1);
%                             uistack(p(iE),'bottom')
                        end
                    end
                else    % e.g. if showEvents is LickIn - BH
                    et = EventTimes(iE,ind);
                    et2 = [];
                    nt2 = [];
                    for iT = 1:NumTrials(iCOMP)
                        tmet = et{iT};
                        tmet(tmet<time(1)) = NaN;
                        tmet(tmet>time(end)) = NaN;
                        et2 = [et2 tmet'];
                        nt2 = [nt2 ones(1,length(tmet))*nt(iT)];
                    end
                    p(iE) = patch([et2-det; et2-det; et2+det; et2+det],[nt2+dnt; nt2-dnt; nt2-dnt; nt2+dnt],'k');
                    try
                        set(p(iE),'FaceColor',par.ShowEventsColors{iE},'EdgeColor','none');
                        uistack(p(iE),'bottom')
                    catch
                        set(p(iE),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.1);
                        uistack(p(iE),'bottom')
                    end
                end
            else    % deals with lick-aligned raster - BH
                et = cell2mat(EventTimes{iE}(ind));
                et(et<time(1))   = NaN;
                et(et>time(end)) = NaN;
                p(iE) = patch([et'-det; et'-det; et'+det; et'+det; et'-det],[nt+dnt; nt-dnt; nt-dnt; nt+dnt; nt+dnt],'k');
%                 l=line([et'; et'],[(nt-0.5); (nt+0.5)],'Color',par.ShowEventsColors{iE},'Linewidth',par.SpikeWidth*16);
                try
                    set(p(iE),'FaceColor',par.ShowEventsColors{iE},'EdgeColor',par.ShowEventsColors{iE});
                    uistack(p(iE),'bottom')
                catch
                    set(p(iE),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.1);
                    uistack(p(iE),'bottom')
                end
            end
        end   % iE
        
        % Axes, labels
        axis([par.window(1) par.window(2) 1-2*eps nt(end)]);
        if (par.RasterPlotNumtrials)
            set(ax(iAX),'YAxisLocation','left','YColor','k','YTick',[1-2*eps nt(end)],'YTickLabel',{'',num2str(nt(end))});
        end
        
        if NumTrials(iCOMP) > 1
            set(ax(iAX),'XTickLabel','','XColor','w');
            set(ax(iAX),'YAxisLocation','left','XColor','w');           
        end
        yl(iAX) = ylabel(TAGS{iCOMP});
        if strcmp(par.PlotZeroLine,'on'),
            line([0 0],ylim(ax(iAX)),'Color',par.PlotZeroLineColor)
        end
        iAX = iAX + 1;
    end   % NumTrials > 0
end   % iCOMP
NewIndex = find(NumTrials~=0);
NumTrials = NumTrials(NewIndex);
NumParts = length(NumTrials);

% Put up panel for PSTH
if strcmpi(par.PSTH,'on')
    fhandle0 = ax;
else
    set(ax(end),'XTickLabelMode','auto','XColor','k','FontName','Ariel','FontSize',12);
end
set(ax,'TickDir','out');
set(yl,'FontName','Ariel','FontSize',10);

% Rearrange
if ~isnan(par.NumTrials2Plot)
    
    % Rearrange equally
    pos = get(ax(1),'Position');
    for i = 2:NumPanels
        set(ax(i),'Position',pos-[0 (i-1)*(pos(4)*(1+par.RasterPanelDistance)) 0 0]);
    end
else   
    posS = get(ax(1),'Pos');
    posE = get(ax(NumParts),'Pos');
    
    % Distance between panels
    deltaPD = ceil(posS(4)*par.RasterPanelDistance*1000) / 1000;
    deltaY = 0.6;
    TotalTrials = sum(NumTrials);
    pbottom(NumParts) = 0.06 + 0.2 + deltaPD;
    for i = NumParts:-1:1
        pheight(i) = max(ceil(deltaY*NumTrials(i)/TotalTrials*100)/100,deltaPD);
        set(ax(i),'Position',[posS(1) pbottom(i) posS(3) pheight(i)]);
        if i > 1
            pbottom(i-1) = pbottom(i) + pheight(i) + deltaPD;
        end
    end   % i
end   % NumTrial2Plot == NaN

% To remove the labels
if par.ShowRasterTAGS == 0
    delete(yl)
end

% Set PSTH Axes
if strcmpi(par.PSTH,'on')
    posE = get(ax(NumParts),'Pos');
    posPSTH = get(ax(end),'Pos');
    if NumPsthPlots == 1
        try
            set(ax(end),'Pos',[posPSTH(1) posPSTH(2)+par.PsthPanelDistance posPSTH(3) posE(2)-(posPSTH(2)+par.PsthPanelDistance)-0.01])
        catch
            warning('plot_raster2a:setAxis','Error setting PSTH axes.')
        end
    end
end

% Raster ID bars
DXC = par.RasterIDbarDX*2;
if par.RasterIDbarWidth ~= 0
    if (par.RasterPlotNumtrials)
        set(ax(1:NumParts),'YAxisLocation','right');
    else    
        set(ax(1:NumParts),'YColor','w');
    end
    set(yl,'String','');
    axb = zeros(1,NumParts);
    pb = zeros(1,NumParts);
    for iP = 1:NumParts
        pos = get(ax(iP),'Position');
        axb(iP) = axes('Position',[pos(1)-par.RasterIDbarDX-DXC pos(2) DXC pos(4)]);
        jP = NewIndex(iP);
        pb(iP) = patch([0 0 0.5 0.5],[1 0 0 1],par.Colors{jP});
        if isfield(par,'Colors2')
            set(pb(iP),'EdgeColor',par.Colors2{jP},'LineWidth',2);
        else
            set(pb(iP),'EdgeColor','none');
        end
        axis fill
        axis off
    end
end

% -------------------------------------------------------------------------
function child_ax = splitaxes(fhandle0,numchilds)

% Splits an axis into multiple axes.
parent_ax = get(fhandle0,'Position');
width = parent_ax(:,4) / numchilds;
newaxpos = repmat(parent_ax,numchilds,1);
newaxpos(:,4) = width;
newaxpos(:,2) = fliplr(newaxpos(1,2)+[0:length(newaxpos(:,2))-1]*width);
child_ax = zeros(1,numchilds);
for iA = 1:numchilds,
    child_ax(iA) = subplot('Position',newaxpos(iA,:));
end