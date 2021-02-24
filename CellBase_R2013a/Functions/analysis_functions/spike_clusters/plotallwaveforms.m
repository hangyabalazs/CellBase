function [weds, mean_all, H_all] = plotallwaveforms(cellid,varargin)
%PLOTALLWAVEFORMS   Plot waveforms.
%   [W, M] = PLOTALLWAVEFORMS(CELLID) plots individual and average
%   waveforms for all tetrode channels. If a 'MAXNUM', VALUE optional input
%   parameter-value pair is provided, then only a limited number of
%   randomly chosen spikes are plotted. All (or selected) waveforms are
%   returned in W; average waveforms are returned in M.
%
%   See also PLOTWAVEFORMS.

%   Balazs Hangya
%   hangya.balazs@koki.hu
%   15-Dec-2020

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'maxnum',Inf,@isnumeric)   % maximal number of spikes to plot
parse(prs,cellid,varargin{:})
g = prs.Results;

% All spikes
wave_all = extractSpikeWaveforms(cellid,'all','chans','all');  % get all waveforms

% Downsample
if g.maxnum < Inf
    nm_evoked = size(wave_all,1);
    rp = randperm(min(g.maxnum,nm_evoked));
    weds = wave_all(rp,:,:);   % downsample for plotting
else
    weds = wave_all;
end

% Average waveforms
mean_all = squeeze(nanmean(wave_all,1));

% Plot waveforms
warning off
H_all = figure('Position',[624 126 1092 852]);
H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(weds,3)]);
for sp = 1:4
    hold(H(sp),'on')
    plot(H(sp),transpose(squeeze(weds(:,sp,:))),'Color',[0.7 0.7 0.7])
    plot(H(sp),transpose(mean_all(sp,:)),'Color','k','LineWidth',6)
end
warning backtrace
title('Spike shape')