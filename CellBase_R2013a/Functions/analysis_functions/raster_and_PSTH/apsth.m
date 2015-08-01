function [psth psth_sd] = apsth(binraster,dt)
%APSTH   Spike density function with adaptive Gaussian kernel.
%	[PSTH_ACONV PSTH_SD] = APSTH(BINRASTER,DT) computes the spike density
%	function with an adaptive Gaussian kernel optimized for local
%	probability of spiking. For near-zero spiking probability, the kernel
%	is flat resulting in a moving average. For a probability of 1, the
%	kernel converges to a Dirac delta. Note that the mapping of
%	probabilities to kernel width is arbitrary between the extremes. In
%	this implementation, the ALPHA parameter of the GAUSSWIN function
%	(which is inversly related to the SD of the Gaussian window) is
%	linearly mapped to probabilities.
%
%   See also BINRASTER2APSTH and GAUSSWIN.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   25-May-2011

%   Edit log: BH 5/25/11, 5/18/13

% Trial number and epoch length
[tno tl] = size(binraster);

% Merged spike train
[x0 allspks] = find(binraster==1);
ts = sort(allspks)';

% % Calculate adaptive SDF with variable Gaussian Kernel #1 (slightly faster 
% % but no convolved raster)
% spno = length(ts);
% agvd = zeros(1,tl);
% tno_true = nansum(~isnan(binraster));
% prob = nanmean(binraster) / (dt * 1000);
% for t = 1:spno
%     spi = ts(t);
%     tspt = zeros(1,tl);
%     tspt(spi) = 1;
%     if prob(spi) > 1
%         keyboard
%     end
%     wbh = gausswin(9,prob(spi)*50);   % kernel
%     wbh = wbh / sum(wbh);
%     agvd = agvd + filtfilt(wbh,1,tspt);   % convolution from both directions
% end
% psth_aconv = [agvd./tno_true] / dt * 1000;   % SDF
% 
% H = figure;  % Plot SDF
% plot(time,psth_aconv,'k')
% xlim([time(1) time(end)])

% Calculate adaptive SDF with variable Gaussian Kernel #2
agvd = zeros(tno,tl);
prob = nanmean(binraster) / (dt * 1000);
mltp = binraster;
mltp(~isnan(mltp)) = 1;
for k = 1:tno   % convolve trial-wise
    spks = find(binraster(k,:)==1);
    spno = length(spks);
    for t = 1:spno
        spi = spks(t);
        tbinraster = zeros(1,tl);
        tbinraster(spi) = 1;
        if prob(spi) > 1 / (dt * 1000)
            keyboard
        end
        wbh = gausswin(9,prob(spi)*50);
        wbh = wbh / sum(wbh);
        agvd(k,:) = agvd(k,:) + filtfilt(wbh,1,tbinraster);
    end
end
psth = nanmean(agvd.*mltp) / dt;
psth_sd = nanstd(agvd.*mltp/dt);
psth_err = nanstd(agvd.*mltp/dt) / sqrt(tno);