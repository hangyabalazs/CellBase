function mclust_projection_preview(sessionID,varargin)
%MCLUST_PROJECTION_PREVIEW   Plot feature data for clusters.
%   MCLUST_PROJECTION_PREVIEW(SESSIONID) plots feature data for a recording
%   session (SESSIONID). Spikes are plotted in Energy feature space and
%   light-evoked spikes are overlayed in blue. Time period for light-evoked
%   activity is either selected automatically (see FINDSTIMPERIOD) or fixed
%   from 1 to 7 ms post stimulus (default). The default behaviors can be
%   modified using the following optional input arguments (parameter, value
%   pairs, with default values):
%       'stim_period', [0.001 0.007] - start and end of stimulation period
%       after each
%           light pulse; set to empty for auto-detection
%       'feature_names', 'Energy' - features for which feature data are
%           plotted; character or cell array
%       'marker', '+' - marker for the scatter plots 'marker_size', 2 -
%       marker size for the scatter plots 'usefastplot', true - use fast
%       plotting method (downsample points
%           to plot only one point per pixel; appears the same); faster,
%           but zoom is not implemented (for saving in image formats, e.g.
%           pdf or jpg); if false, full data is plotted (for viewing or
%           saving fig format)
%       'stimonly', true - only spikes from the beginning of the first to
%           the end of the last stimulation protocol are selected for
%           plotting; if false, all spikes are included
%       'plotlightspikes', true - if true, light-evoked spikes are
%           superimposed
%
%   Example:
%   mclust_projection_preview({'n084' '150622a'},'feature_names',...
%       {'Amplitude','Energy'},'stim_period',[0.002 0.004]);
%
%   See also PLOT_MCLUST_PROJECTIONS2.

%   Balazs Hangya
%   Institute of Experimental Medicine, Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu
%   18-Oct-2016

% Input arguments
prs = inputParser;
addRequired(prs,'sessionID',@issessionid)
addParamValue(prs,'stim_period',[0.001 0.007],@isnumeric)   % start and end point for the interval of light-evoked spikes
addParamValue(prs,'feature_names',{'Energy'},@(s)iscell(s)|ischar(s))  % names of MClust features to use
addParamValue(prs,'marker','+')  % marker for the plots
addParamValue(prs,'marker_size',2,@isnumeric)  % marker size for the plots
addParamValue(prs,'usefastplot',true,@(s)islogical(s)|ismember(s,[0 1]))  % fast plotting (no zoom)
addParamValue(prs,'plotlightspikes',true,@(s)islogical(s)|ismember(s,[0 1]))  % plot light-evoked spikes
parse(prs,sessionID,varargin{:})
g = prs.Results;
if ischar(g.feature_names)
    g.feature_names = {g.feature_names};
end

% Load trial events (unsynchronized)
cbdir = getpref('cellbase','datapath');
sessionpath = fullfile(cbdir,sessionID{1,1},sessionID{1,2});
[EventTimeStamps, EventIDs, Nttls, Extras, EventStrings NlxHeader] = ...
    Nlx2MatEV([sessionpath '\Events.nev'],[1 1 1 1 1],1,1,1);

% Load spikes from Ntt file
dr = dir(sessionpath);
tt = getpref('cellbase','cell_pattern');
% clsinx = strncmp({dr.name},tt,2);
clsinx = cellfun(@(s)~isempty(s),regexp({dr.name},[tt '\w*.ntt']));
cls = {dr(clsinx).name};
numTetrodes = length(cls);   % number of tetrodes

for iT = 1:numTetrodes   % tetrode loop
    Nttfn = fullfile(sessionpath,cls{iT});   % spike data file name
    
    all_spikes = LoadTT_NeuralynxNT(Nttfn);
    TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
    all_spikes = all_spikes * TIMEFACTOR;
    numSpikes = length(all_spikes);
    
    % Load feature data for tetrode.
    numFeatures = length(g.feature_names);   % number of spikeshape features
    FeatureData = nan(numFeatures,numSpikes,4);
    for k = 1:numFeatures   % feature loop
        basename = [tt num2str(iT)];
        propfn = [basename '_' g.feature_names{k}];   % name of feature file (e.g. TT1_Amplitude)
        propfn_path = [sessionpath filesep 'FD'];   % where the feature file can be found
        if ~isdir(propfn_path)
            propfn_path = sessionpath;
        end
        propfn_full = [propfn_path filesep propfn];   % full path of feature file
        try
            wf_prop = load([propfn_full '.fd'],'-mat');     % load feature file
        catch    %#ok<CTCH>
            disp('Calculating missing feature data.')       % calculate feature file if it was not found
            valid_channels = [1 1 1 1];   % calculate features for all channels
            calculate_features(sessionpath,propfn_path,g.feature_names(k),basename,valid_channels)
            wf_prop = load([propfn_full '.fd'],'-mat');     % load feature file
        end
        FeatureData(k,:,:) = wf_prop.FeatureData;
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
        inx = Nttls==16384 | Nttls==16;
        light_ts = EventTimeStamps(inx) / 10^2 * TIMEFACTOR;
        f_ts = all_spikes;
        show_ts = [];
        for k = 1:length(light_ts)
            show_ts = [show_ts; find(f_ts>light_ts(k)+lim1&f_ts<light_ts(k)+lim2)];
        end
        LightSpikes = show_ts;   % indices for light spikes
    end
    
    % Plot
    pcmb = allcomb(1:size(FeatureData,1),1:size(FeatureData,3));
    cmb = flipud(combnk(pcmb(:,1)*10+pcmb(:,2),2));
    NumComb = size(cmb,1);   % number of projections
    for k = 1:NumComb   % projection loop
        fst = [floor(cmb(k,1)/10) mod(cmb(k,1),10)];
        scnd = [floor(cmb(k,2)/10) mod(cmb(k,2),10)];
        xdata = squeeze(FeatureData(fst(1),:,fst(2)));
        ydata = squeeze(FeatureData(scnd(1),:,scnd(2)));
        
        % Open figure
        namestr = [g.feature_names{fst(1)} num2str(fst(2)) '_'...
            g.feature_names{scnd(1)} num2str(scnd(2))];
        HS.(namestr) = figure;  % 'Position',[624 126 1092 852]
        hold on
        axis([min(xdata) max(xdata)+1 min(ydata) max(ydata)+1])
        xlabel([g.feature_names{fst(1)} ': ' num2str(fst(2))])
        ylabel([g.feature_names{scnd(1)} ': ' num2str(scnd(2))])
        
        % Plot all spikes
        if g.usefastplot
            fastplot(xdata,ydata,[0 0 0],g.marker,g.marker_size);
        else
            slowplot(xdata,ydata,[0 0 0],g.marker,g.marker_size);
        end
        
        % Plot light-evoked spikes
        if g.plotlightspikes
            xdatai = squeeze(FeatureData(fst(1),LightSpikes,fst(2)));
            ydatai = squeeze(FeatureData(scnd(1),LightSpikes,scnd(2)));
            if g.usefastplot
                fastplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size+5);
            else
                slowplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size+5);
            end
        end
    end
    
    % Print figures to pdf
    pdfname = fullfile(sessionpath,['CLUSTERPREVIEW_' sessionID{1} '_' ...
        sessionID{2} '_TT' num2str(iT) '.tiff']);
    writefigs(HS,pdfname)
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

% -------------------------------------------------------------------------
function calculate_features(sessionpath,propfn_path,feature_names,basename,valid_channels)

% Create MClust variables
global MClust_FDdn
global MClust_ChannelValidity
global MClust_NeuralLoadingFunction
global MClust_TText
global MClust_FDext
global MClust_TTdn
global MClust_TTfn

MClust_FDext = '.fd';
MClust_TText = '.ntt';
MClust_TTfn = basename;
[t1, t2] = strtok(fliplr(which('MClust')),filesep);
MClust_Directory = fliplr(t2);
MClust_FDdn = propfn_path;
MClust_TTdn = sessionpath;
MClust_ChannelValidity = valid_channels;
MClust_NeuralLoadingFunction = char([MClust_Directory 'LoadingEngines\LoadTT_NeuralynxNT']);

% Calculate features
CalculateFeatures(basename,feature_names)

% -------------------------------------------------------------------------
function writefigs(H,pdfname)
% Write figures to pdf
 
% Append to pdf
if isstruct(H)  % H is either a struct with figure handles, or a single fig. handle
    fls = fieldnames(H);
    for fs = 1:length(fls)
        h = H.(fls{fs});
        if ishandle(h)
            export_fig(h,'-append',pdfname, '-painters');
            close(h)
%             fh = ['-f' num2str(h.Number)];
%             print(fh,'-dpdf',pdfname)  % write to pdf
        end
    end
else
    export_fig(H,'-append',pdfname, '-opengl');  % write to pdf
    close(H)
%     fh = ['-f' num2str(H.Number)];
%     print(fh,'-dpdf',pdfname)  % write to pdf
end