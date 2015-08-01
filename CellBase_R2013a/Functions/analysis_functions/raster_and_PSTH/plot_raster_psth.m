function HS = plot_raster_psth(cellid,varargin)
%PLOT_RASTER_PSTH   Plot raster plots and PSTHs.
%   HS = PLOT_RASTER_PSTH(CELLID) plots raster plots and PSTHs aligend to
%   'BurstOn' and/or 'PulseOn' events. Figure handles are returned in HS.
%   Which rasters/PSTHs are plotted is controlled by the following optional
%   input argument parameter-value pairs (with default values):
%       'BurstOn', true - if true, plot 'BurstOn' raster and PSTH
%       'PulseOn', false - if true, plot 'PulseOn' raster and PSTH
%
%   Output (fields of HS struct):
%       H_BurstOn - figure handle for 'BurstOn' raster and PSTH
%       H_PulseOn - figure handle for 'PulseOn' raster and PSTH
%
%   See also VIEWCELL2B and PLOT_RASTER2A.

%   Edit log: BH 5/9/12

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'BurstOn',true,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying 'BurstOn' rasters and PSTHs
addParamValue(prs,'PulseOn',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying 'PulseOn' rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;

% Set input parameters for 'viewcell2b'
SEvent = 'BurstOff';
win = [-0.2 0.5];
partsPO = 'all';
partsBO = '#BurstNPulse';
dt = 0.001;
sigma = 0.001;
PSTHstd = 'on';
ShEvent = {{'PulseOn','PulseOff','BurstOff'}};
ShEvColors = hsv(length(ShEvent{1}));
ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);

% Plot raster plot and PSTH for 'BurstOn'
if g.BurstOn
    HS.H_BurstOn = figure('Position',[97 163 1533 815]);
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName','BurstOn','SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',gcf,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',partsBO,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% Plot raster plot and PSTH for 'PulseOn'
if g.PulseOn
    HS.H_PulseOn = figure('Position',[97 163 1533 815]);
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName','PulseOn','SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',gcf,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',partsPO,...
        'EventMarkerWidth',0,'PlotZeroLine','off','Num2Plot',1000)
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end