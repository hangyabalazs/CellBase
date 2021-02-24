function f0 = plot_timecourse(time,PSTH,PSTHse,params,varargin)
%PLOT_TIMECOURSE   Plot PSTH.
%   F0 = PLOT_TIMECOURSE(TIME,PSTH,PSTHSE) plots PSTH as a function of TIME
%   and uses PSTHSE to generate error shade around the line; the figure
%   handle F0 is returned. Additional parameters can be passed as
%   parameter, value pairs (with default values):
%       'TitleStr', '' - title string
%       'XLabel', 'Time' - x axis label
%       'YLabel', 'Response' - y axis label
%       'ShowPsthYLabel', 1 - controls y axis labeling
%       'Legend', '' - figure legend
%       'PSTHstd', 'off' - controls error shading
%       'PSTHlinewidth', 2 - line width
%       'LineStyle', '' - line style
%       'PlotZeroLine', 'on' - controls whether to plot a line at zero
%       'PlotZeroLineColor', [0 0 0] - color for the line at zero.
%
%   See also VIEWCELL2B and PLOT_RASTER2A.

%   Edit log: BH 4/15/13

% Default arguments
default_args = params;
default_args.TitleStr = '';
default_args.XLabel = 'Time';
default_args.YLabel = 'Response';
default_args.ShowPsthYLabel = 1;
default_args.Legend = '';
default_args.PSTHstd = 'off';
default_args.PSTHlinewidth = 2;
default_args.LineStyle = '';
default_args.PlotDashedTime = NaN;
default_args.DashedLineStyle = '--';
default_args.ClearFig = 'off';
default_args.PlotZeroLine = 'on';
default_args.PlotZeroLineColor = [0 0 0];
[par,error] = parse_args(default_args,params);  % overwrite defaults with passed arg structure
[par,error] = parse_args(par,varargin{:});      % overwrite all with passed individual args

% Restrict to the given time frame
pos2disp = restrict2(time,par.window(1),par.window(2));

% Calculate extrema
MX = max(max(PSTH(:,pos2disp)));
MN = min(min(PSTH(:,pos2disp)));
RANGE = (MX-MN);
MX = MX + RANGE * 0.25;
MN = MN - RANGE * 0.25;
alim = [par.window(1) par.window(2) MN MX+200*eps];

% Focus figure 
if strcmp(get(par.FigureNum,'Type'),'axes'),
    f0 = focusfigure(par.FigureNum);
    cla;
else
    f0 = focusfigure(par.FigureNum,'create');
end
if strcmpi(par.ClearFig,'on')
    clf;
end
hold on;

% Plot
NumPlots = [];
for i = 1:size(PSTH,1)
    if strcmpi(par.PSTHstd,'on')
        if ~isempty(par.LineStyle)
            errorshade(time(pos2disp),PSTH(i,pos2disp),PSTHse(i,pos2disp),'LineColor',par.Colors{i},...
                'LineStyle',par.LineStyle{i},'LineWidth',par.PSTHlinewidth,'PlotDashedTime',par.PlotDashedTime,'DashedLineStyle',par.DashedLineStyle);
        else
            errorshade(time(pos2disp),PSTH(i,pos2disp),PSTHse(i,pos2disp),'LineColor',par.Colors{i},'LineWidth',par.PSTHlinewidth,'PlotDashedTime',par.PlotDashedTime); 
        end
    else
        if ~isnan(par.PlotDashedTime)
            posDASHED = nearest(time(pos2disp),par.PlotDashedTime);
            p0(1) = plot(time(pos2disp(1:posDASHED)),PSTH(i,pos2disp(1:posDASHED)));
            p0(2) = plot(time(pos2disp(posDASHED:end)),PSTH(i,pos2disp(posDASHED:end)),par.DashedLineStyle);
        else
            p0 = plot(time(pos2disp),PSTH(i,pos2disp));
        end
        if ~isempty(par.LineStyle)
            set(p0(1),'LineStyle',par.LineStyle{i});
        end
        set(p0,'Color',par.Colors{i},'LineWidth',par.PSTHlinewidth);
    end
    if min(1,sum(~isnan(PSTH(i,:))))
        NumPlots = [NumPlots i];
    end
end

% Labels
l(1) = ylabel(par.YLabel);
l(2) = xlabel(par.XLabel);
li = 2;
if ~isempty(par.TitleStr)
    li = li+1;
    l(li) = title(par.TitleStr);
end
if par.ShowPsthYLabel == 0
    delete(l(1));
end

% Legend
if ~isempty(par.Legend)
    l_leg = legend(par.Legend(NumPlots),'Location','NorthEast');
    set(l_leg,'box','off','FontSize',8,'Color','none');
end
set(gca,'XTickLabelMode','auto','XColor','k','YTickLabelMode','auto','YColor','k','YTickMode','auto');
axis(alim);
setmyplot(gca,l);