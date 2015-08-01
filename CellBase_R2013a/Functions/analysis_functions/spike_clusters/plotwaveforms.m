function out = plotwaveforms(cellid,varargin)
%PLOTWAVEFORMS   Plot spontaneous and/or light-evoked spike shape.
%   PLOTWAVEFORMS(CELLID) plots the waveform of spontaneous and putative
%   light-evoked spikes. It finds the latency of putative light-evoked
%   spikes as peaks of firing rate after light pulses. The period for
%   putative light-evoked spikes is defined by half-hight crossings around
%   the peak (see FINDSTIMPERIOD for details). If there is no detectable
%   activation, 1 and 6 ms are used, which include typical latencies of
%   light-evoked spikes. Spontaneous spikes are extracted from 2 s periods
%   before light bursts. Only 2000 randomly chosen spikes are plotted. This
%   number can be changed using the 'MAXNUM', VALUE optional input
%   parameter-value pair.
%
%   Two optional input parameter value pairs control whether to plot both
%   spontaneous and evoked waveforms or only one of those:
%   HS = PLOTWAVEFORMS(CELLID,'SPONT',TRUE,'EVOKED',TRUE) plots both
%   waveforms and a third plot to compare them. Setting either of the input
%   values to 'false' will prevent plotting that type. Handles of the
%   resulting figures are returned in the output structure HS with the
%   following fieldnames: HS.H_spont for spontaneous waveforms, HS.H_evoked
%   for light-evoked waveforms and HS.H_compare for the compound figure.
%
%   HS = PLOTWAVEFORMS(CELLID,'CORRELATION','TRUE') also returns waveform
%   correlation for the largest channel (see SPIKESHAPECORR) in HS.R field
%   of the output structure.
%
%   Optional input parameter-value paris with default values:
%       'spont', true - plot spontaneous waveform
%       'evoked', true - plot light-evoked waveform
%       'maxnum', 2000 - maximum number of spikes to plot
%       'correlation', false - return waveform correlation
%       'stim_period', [] - start and end of stimulation period after each
%           light pulse
%
%   Output (fields of HS struct):
%       H_spont - figure handle for spontaneous spike waveform
%       H_evoked - figure handle for evoked spike waveform
%       H_compare - figure handle for comparison plot
%       R - spike shape correlation of spont. and evoked waveform
%       activation_start - start of the detected stimulation period
%       activation_end - end of the detected stimulation period
%
%   See also FINDSTIMPERIOD, ABS2RELTIMES, EXTRACTSEGSPIKES,
%   EXTRACTSPIKEWAVEFORMS and SPIKESHAPECORR.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   07-May-2012

%   Edit log: BH 5/7/12

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'spont',true,@(s)islogical(s)|ismember(s,[0 1])) % plot spontaneous waveform
addParamValue(prs,'evoked',true,@(s)islogical(s)|ismember(s,[0 1])) % plot light-evoked waveform
addParamValue(prs,'maxnum',2000,@isnumeric)   % maximal number of spikes to plot
addParamValue(prs,'correlation',false,@(s)islogical(s)|ismember(s,[0 1])) % calculate waveform correlation
addParamValue(prs,'stim_period',[],@isnumeric)   % start and end point for the interval of light-evoked spikes
parse(prs,cellid,varargin{:})
g = prs.Results;

% Evoked spikes
if g.evoked || g.correlation
    ST = loadcb(cellid,'STIMSPIKES');  % load stimulus spikes (prealigned)
    if isequal(findcellstr(ST.events(:,1),'PulseOnS'),0)  % if prealignment to stim. period has not been saved yet
        if isempty(g.stim_period)   % latency of stimulated spikes
            [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
            if isnan(lim1) || isnan(lim2)
                lim1 = 0.001;   % no activation detected
                lim2 = 0.006;
            end
        else  % stimulation period passed as input
            lim1 = g.stim_period(1);
            lim2 = g.stim_period(2);
        end
        tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'stim','PulseOn');   % convert period to epochs relative to pulses
        selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find putative stimualated spikes
    else  % if prealignment to stim. period has already been saved
        trigger_pos = findcellstr(ST.events(:,1),'PulseOnS');
        pse = ST.event_stimes{trigger_pos};
        selts_evoked = rel2abstimes(cellid,pse,'stim','PulseOn');  % convert to absolute spike times
        lim1 = ST.events{trigger_pos,4}(1);  % boundaries of stim. period used for prealigning spikes
        lim2 = ST.events{trigger_pos,4}(2);
    end
    wave_evoked = extractSpikeWaveforms(cellid,selts_evoked,'chans','all');  % get waveforms for the extracted spikes

    % Downsample
    nm_evoked = size(wave_evoked,1);
    rp = randperm(min(g.maxnum,nm_evoked));
    weds = wave_evoked(rp,:,:);   % downsample for plotting
    
    % Average waveforms
    mean_evoked = squeeze(nanmean(wave_evoked,1));
end

% Spontaneous spikes
if g.spont || g.correlation
    tsegs_spont = rel2abstimes(cellid,[-2,0],'stim','BurstOn');   % extract 2s periods before bursts
    selts_spont = extractSegSpikes(cellid,tsegs_spont);     % extract spontaneous spikes
    wave_spont = extractSpikeWaveforms(cellid,selts_spont,'chans','all');    % get waveforms for the extracted spikes
    
    % Downsample
    nm_spont = size(wave_spont,1);
    rp = randperm(min(g.maxnum,nm_spont));
    wsds = wave_spont(rp,:,:);
    
    % Average waveforms
    mean_spont = squeeze(nanmean(wave_spont,1));
end

% Plot light-evoked waveforms
if g.evoked
    out.H_evoked = figure('Position',[624 126 1092 852]);
    H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
    for sp = 1:4
        hold(H(sp),'on')
        plot(H(sp),transpose(squeeze(weds(:,sp,:))))
    end
    title('Light-evoked spike shape')
end

% Plot spontaneous waveforms
if g.spont
    out.H_spont = figure('Position',[624 126 1092 852]);
    H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
    for sp = 1:4
        hold(H(sp),'on')
        plot(H(sp),transpose(squeeze(wsds(:,sp,:))))
    end
    title('Spontaneous spike shape')
end

% Compare waveforms
if g.evoked && g.spont
    out.H_compare = figure('Position',[624 126 1092 852]);
    H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
    for sp = 1:4
        hold(H(sp),'on')
        plot(H(sp),transpose(squeeze(wsds(:,sp,:))),'Color',[0.9 0.9 0.9])
        plot(H(sp),transpose(mean_spont(sp,:)),'Color','k','LineWidth',6)
        plot(H(sp),transpose(mean_evoked(sp,:)),'Color',[0 153 255]/255,'LineWidth',2)
    end
    title('Compare spont. and light-evoked spike shape')
end

% Spike shape correlation
if g.correlation
    if length(selts_evoked) > 1 && length(selts_spont) > 1   % if there are waveforms to calculate corr. for
        mx = maxchannel(wave_spont);     % mx: largest channel
        mnmx_spont = nanmean(squeeze(wave_spont(:,mx,:)));
        mnmx_evoked = nanmean(squeeze(wave_evoked(:,mx,:)));
        sr = 32552;     % DigiLynx sampling rate
        rng = round(0.00075*sr);    % number of data points in 750 us (default censored period of DigiLynx)
        pr = corrcoef(mnmx_spont(1:rng),mnmx_evoked(1:rng));
        out.R = pr(1,2);
    else
        out.R = NaN;
    end
end

% Additional output
if exist('lim1','var')
    out.activation_start = lim1;
    out.activation_end = lim2;
end

% -------------------------------------------------------------------------
function mx = maxchannel(wv)

% Find largest channel
mean_wv = squeeze(nanmean(wv,1));   % mean waveform
amx = max(max(mean_wv));     % absolut maximum of mean waveforms
[mx my] = find(mean_wv==amx,1,'first');     %#ok<NASGU> % mx: largest channel