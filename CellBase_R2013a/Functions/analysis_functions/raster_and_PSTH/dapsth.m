function [psth psth_sd] = dapsth(binraster,dt)
%DAPSTH   Spike density function with 'doubly adaptive' Gaussian kernel.
%	[PSTH_ACONV PSTH_SD] = DAPSTH(BINRASTER,DT) computes the spike density
%	function with an adaptive Gaussian kernel optimized for local
%	probability of spiking. For near-zero spiking probability, the kernel
%	is flat resulting in a moving average. For a probability of 1, the
%	kernel converges to a Dirac delta. Note that the mapping of
%	probabilities to kernel width is arbitrary between the extremes. In
%	this implementation, the SD of the Gaussian probability density
%	function is linearly mapped on the reciprocal of the probabilities.
%	Also, length of the integration window is linearly mapped on
%	probabilities, with long window corresponding to low probability.
%
%   See also BINRASTER2DAPSTH and GAUSSWIN.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   05-Jul-2012

%   Edit log: BH 7/5/12, 5/18/13

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
% figure
% A = axes;
% hold on
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
        wnlm = (1 - prob(spi)).^5 * 26 + 4;   % map half-window size from +-26 to +-4 (prob. 0 to prob. 1); 
            % exponent sharpens the transition
        wnlm = max(wnlm,1);   % prevent empty window for dt < 0.001
        shf = 0.5;   % 'sharpness factor': scaler that sharpens the kernel for higher probabilities (0.5 - sharp)
        gsd = shf / prob(spi) - 1;   % SD of the Gaussian window
        gsd = max(gsd,0.1);   % normpdf function numerically fails for smaller SD values; with 0.1 it's already practically Dirac delta 
        gsd = min(gsd,10000);   % normpdf function numerically fails for higher SD values; with 10000 it's already practically moving average 
        wbh = normpdf(-wnlm:wnlm,0,gsd);   % map Gaussian windows on probabilities from moving average 
            % (prob. 0) to Dirac delta (prob. 1)
        wbh = wbh / sum(wbh);
%         plot(wbh,'Color',rand(1,3))
        agvd(k,:) = agvd(k,:) + filtfilt(wbh,1,tbinraster);
    end
end
psth = nanmean(agvd.*mltp) / dt;
psth_sd = nanstd(agvd.*mltp/dt);
psth_err = nanstd(agvd.*mltp/dt) / sqrt(tno);

% % Plot
% H1 = figure;
% subplot(211)
% imagesc(agvd);
% % set(gca,'XTick',[])
% subplot(212)
% errorshade(1:length(psth),psth,psth_err,'LineColor','k','LineWidth',2,'ShadeColor','grey')
% % xlim([time(1) time(end)])   
% axes('Position',[0.7 0.35 0.2 0.2])
% Ls = findobj(allchild(A),'type','line');   % copyobj runs to Matlab bug
% for k = 1:length(Ls)
%     plot(get(Ls(k),'XData'),get(Ls(k),'YData'),'Color',get(Ls(k),'Color'))
%     hold on
% end