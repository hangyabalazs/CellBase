function HS = plot_mclust_projections2(cellid,varargin)
%PLOT_MCLUST_PROJECTIONS2   Plot feature data for clusters.
%   PLOT_MCLUST_PROJECTIONS2(CELLID) plots feature data for a tagged cell
%   (CELLID) and it's tetrode pairs in all projections. Tagged cluster is
%   shown in orange and light-evoked spikes are overlayed in blue. Time
%   period for light-evoked activity is selected automatically (see
%   FINDSTIMPERIOD). Only spikes from the beginning of the first to the end
%   of the last stimulation protocol are included. These default behaviors
%   can be modified using the following optional input arguments
%   (parameter, value pairs, with default values):
%       'stim_period', [] - start and end of stimulation period after each
%           light pulse
%       'feature_names', 'Energy' - features for which feature data are
%           plotted; character or cell array
%       'marker', '+' - marker for the scatter plots
%       'marker_size', 2 - marker size for the scatter plots
%       'usefastplot', true - use fast plotting method (downsample points
%           to plot only one point per pixel; appears the same); faster,
%           but zoom is not implemented (for saving in image formats, e.g.
%           pdf or jpg); if false, full data is plotted (for viewing or
%           saving fig format)
%       'MClust_PlotSkipPoints', 0 - disregard the largest N values on x 
%           and y axes for plotting (prevent rescaling by large noise)
%       'stimonly', true - only spikes from the beginning of the first to
%           the end of the last stimulation protocol are selected for
%           plotting; if false, all spikes are included
%       'plotlightspikes', true - if true, light-evoked spikes are 
%           superimposed
%
%   HS = PLOT_MCLUST_PROJECTIONS2(CELLID,...) returns the handles of the
%   figures in HS struct. The fields of HS are named according to the
%   features and channels: HS.(XFeatureXChannel_YFeatureYChannel), e.g.
%   HS.Amplitude1_Energy4.
%
%   Example:
%   plot_mclust_projections2(cellid,'feature_names',{'Amplitude','Energy'},...
%        'stim_period',[0.002 0.004]);
%
%   See also PLOTWAVEFORMS.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   08-May-2012

%   Edit log: SPR 12/28/11; BH 5/8/12

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'stim_period',[],@isnumeric)   % start and end point for the interval of light-evoked spikes
addParamValue(prs,'feature_names',{'Energy'},@(s)iscell(s)|ischar(s))  % names of MClust features to use
addParamValue(prs,'marker','+')  % marker for the plots
addParamValue(prs,'marker_size',2,@isnumeric)  % marker size for the plots
addParamValue(prs,'usefastplot',true,@(s)islogical(s)|ismember(s,[0 1]))  % fast plotting (no zoom)
addParamValue(prs,'MClust_PlotSkipPoints',0,@isnumeric)  % disregard the largest N values on x and y axes
addParamValue(prs,'stimonly',true,@(s)islogical(s)|ismember(s,[0 1]))  % restrict to period between first and last light pulse
addParamValue(prs,'plotlightspikes',true,@(s)islogical(s)|ismember(s,[0 1]))  % plot light-evoked spikes
parse(prs,cellid,varargin{:})
g = prs.Results;
if ischar(g.feature_names)
    g.feature_names = {g.feature_names};
end

% Load stim events to get Pulse onset times.
ST = loadcb(cellid,'Stimevents');
pon = ST.PulseOn(~isnan(ST.PulseOn));

% Load spikes from Ntt file.
% Nttfn = cellid2fnames(cellid,'Ntt');
Nttfn = cellid2fnames(cellid,'OE');
% all_spikes = LoadTT_NeuralynxNT(Nttfn);
all_spikes = LoadTT_Intan(Nttfn);
TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
all_spikes = all_spikes * TIMEFACTOR;
val_spk_i = [find(all_spikes >= pon(1),1,'first') ...
    find(all_spikes <= pon(end),1,'last')]; % consider spikes only within the stimulation protocol to account for drift
nspk = length(all_spikes);
spk = loadcb(cellid,'Spikes');
[junk,junk2,tagged_cell_inx] = intersect(spk,all_spikes);  %#ok<*ASGLU> % get indices for the cell
if ~isequal(junk,spk)  % check if all files have appropriate time stamps
    error('plot_mclust_projections:SpikeTimeMismatch','Mismatch between saved spike times and Ntt time stamps.')
end
if g.stimonly   % restrict to stimulation epoch
    tagged_cell_inx = tagged_cell_inx(tagged_cell_inx>val_spk_i(1)&tagged_cell_inx<val_spk_i(2));
end
tagged_cell_i = zeros(nspk,1);
tagged_cell_i(tagged_cell_inx) = 1;

% Cells on same tetrode including cellid.
[NumCell,tetpartners] = tetrodepairs(cellid);
tag_cell = strmatch(cellid,tetpartners);

% Spikes from each cell have an index. Spikes from noise have index 0.
cell_i = zeros(nspk,1);
for iCell = 1:NumCell
    spk = loadcb(tetpartners(iCell),'Spikes'); % load spike times.
    [junk,junk2,cell_inx] = intersect(spk,all_spikes); % get indices for the cell
    if g.stimonly   % restrict to stimulation epoch
        cell_inx = cell_inx(cell_inx>val_spk_i(1)&cell_inx<val_spk_i(2));
    end
    cell_i(cell_inx) = iCell;
end

% Load feature data for tetrode.
[r,s,t] = cellid2tags(cellid);
for k = 1:length(g.feature_names)
    prop = [g.feature_names{k} '.fd'];
    propfn = [getpref('cellbase','cell_pattern') num2str(t) '_' prop];
    sessionpath = cellid2fnames(cellid,'sess');
    propfn_path = [sessionpath filesep 'FD'];   % where the feature file can be found
    if ~isdir(propfn_path)
        propfn_path = sessionpath;
    end
    propfn_path = fullfile(propfn_path,propfn);
    wf_prop = load(propfn_path,'-mat');
    FeatureData(k,:,:) = wf_prop.FeatureData; %#ok<AGROW>
end

% Light-evoked spikes
if g.plotlightspikes
    
    % Latency of stimulated spikes
    if isempty(g.stim_period)
        [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
        if isnan(lim1) || isnan(lim2)
            lim1 = 0.001;   % no activation detected
            lim2 = 0.006;
        end
    else
        lim1 = g.stim_period(1);
        lim2 = g.stim_period(2);
    end
    
    % Evoked spikes
    tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'stim','PulseOn');   % convert period to epochs relative to pulses
    selts_evoked = extractSegSpikes(all_spikes,tsegs_evoked);   % find putative stimualated spikes
    [junk,junk2,evoked_cell_inx] = intersect(selts_evoked,all_spikes); % get indices for light-evoked spikes
end

% Plot
cmp = hsv(NumCell) / 4 + 0.75;
pcmb = allcomb(1:size(FeatureData,1),1:size(FeatureData,3));
cmb = flipud(combnk(pcmb(:,1)*10+pcmb(:,2),2));
NumComb = size(cmb,1);
for k = 1:NumComb
    fst = [floor(cmb(k,1)/10) mod(cmb(k,1),10)];
    scnd = [floor(cmb(k,2)/10) mod(cmb(k,2),10)];
    xdata = squeeze(FeatureData(fst(1),cell_i==0,fst(2)));
    ydata = squeeze(FeatureData(scnd(1),cell_i==0,scnd(2)));
    
    % Open figure
    namestr = [g.feature_names{fst(1)} num2str(fst(2)) '_'...
        g.feature_names{scnd(1)} num2str(scnd(2))];
    HS.(namestr) = figure('Position',[624 126 1092 852]);
    hold on
    sx = sort(xdata);
    mnx = sx(1);
    mxx = sx(end-g.MClust_PlotSkipPoints);   % disregard the largest N values on x axis
    sy = sort(ydata);
    mny = sy(1);
    mxy = sy(end-g.MClust_PlotSkipPoints);   % disregard the largest N values on y axis

    axis([mnx mxx+1 mny mxy+1])
    xlabel([g.feature_names{fst(1)} ': ' num2str(fst(2))])
    ylabel([g.feature_names{scnd(1)} ': ' num2str(scnd(2))])
    
    % Plot noise spikes
    if g.usefastplot
        fastplot(xdata,ydata,[0.8 0.8 0.8],g.marker,g.marker_size);
    else
        slowplot(xdata,ydata,[0.8 0.8 0.8],g.marker,g.marker_size);
    end
    
    % Plot all clusters
    for iC = 1:NumCell
        xdatai = squeeze(FeatureData(fst(1),cell_i==iC,fst(2)));
        ydatai = squeeze(FeatureData(scnd(1),cell_i==iC,scnd(2)));
        if g.usefastplot
            fastplot(xdatai,ydatai,cmp(iC,:),g.marker,g.marker_size);
        else
            slowplot(xdatai,ydatai,cmp(iC,:),g.marker,g.marker_size);
        end
    end
    
    % Plot tagged cluster
    xdatai = squeeze(FeatureData(fst(1),tagged_cell_i==1,fst(2)));
    ydatai = squeeze(FeatureData(scnd(1),tagged_cell_i==1,scnd(2)));
    if g.usefastplot
        fastplot(xdatai,ydatai,[255 204 0]/255,g.marker,g.marker_size);
    else
        slowplot(xdatai,ydatai,[255 204 0]/255,g.marker,g.marker_size);
    end
    
    % Plot light-evoked spikes
    if g.plotlightspikes
        xdatai = squeeze(FeatureData(fst(1),evoked_cell_inx,fst(2)));
        ydatai = squeeze(FeatureData(scnd(1),evoked_cell_inx,scnd(2)));
        if g.usefastplot
            fastplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size+5);
        else
            slowplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size+5);
        end
    end
end

% -------------------------------------------------------------------------
function C = allcomb(A,B)

% Convert to columns
A = A(:);
B = B(:);

% Combinations
as = arrayfun(@(k)horzcat(repmat(A(k),length(B),1),B),1:length(A),'UniformOutput',false);
C = cell2mat(as');

% -------------------------------------------------------------------------
function h = fastplot(x,y,clr,mrk,mrks)

% Reduce number of points
old_units = get(gca,'Units');
set(gca,'Units','pixels')
pos = get(gca,'Position');
xpixels = pos(3) + 1;
ypixels = pos(4) + 1;

xl = xlim;
mnx = xl(1);
mxx = xl(2);
yl = ylim;
mny = yl(1);
mxy = yl(2);
x2 = round((x-mnx)/(mxx-mnx)*(xpixels-1)) + 1;
y2 = round((y-mny)/(mxy-mny)*(ypixels-1)) + 1;
u = unique(x2*100000+y2);
y3 = mod(u,100000);
x3 = (u - y3) / 100000;
x4 = (x3 / xpixels) * (mxx - mnx) + mnx;
y4 = (y3 / ypixels) * (mxy - mny) + mny;

% Plot
h = plot(x4,y4,'.','Color',clr,'Marker',mrk,'MarkerSize',mrks);

% Restore axis units property
set(gca,'Unit',old_units)

% -------------------------------------------------------------------------
function h = slowplot(x,y,clr,mrk,mrks)

h = plot(x,y,'.','Color',clr,'Marker',mrk,'MarkerSize',mrks);