function [p_value Idiff] = nbisstim(cellid)
%NBISSTIM   Assessment of optic tagging.
%   [P I] = NBISSTIM(CELLID) calculates information distances and
%   corresponding p values for light tagging for the cell given in CELLID
%   (see CellBase documentation). Briefly, a spike raster is calculated
%   with 1 ms resolution. The raster is devided to 10 ms time bins and
%   spike latency distribution for first spikes is computed within each
%   bin. Pairwise information divergence measures are calculated for the
%   before-light distributions to form a null-hypothesis distribution for
%   distances. The distances of the first after-light distribution (using
%   'BurstOn' events - see CellBase documentation) from all before-light
%   distributions is calculated and the median of these values is tested
%   against the null-hypothesis distribution. Output arguments P and I
%   correspond to a modified version of Jensen-Shannon divergence (see
%   Endres and Schindelin, 2003)
%
%   NBISSTIM takes a different bin raster (typically aligned to 'BurstOn'
%   events) for baseline distribution than the the one for testing against
%   baseline (typically aligned to 'PulseOn' events). Number of test trials
%   used is maximized to 5000.
%
%   Reference:
%   Endres DM, Schindelin JE (2003) A new metric for probability
%   distributions. IEEE Transactions on Information Theory 49:1858-1860.
%
%   See also TAGGING, STIMES2BINRASTER and JSDIV.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com

% Input argument check
if nargin < 2
    win = [-0.6 0.6];  % time window for bin raster (for Hyun: 1.2)
    dt = 0.001;   % resolution of bin raster in s (for Hyun: 0.0005)
    dsply = 0;   % repress display
end

% Set parameters and load CellBase variables
EventName1 = 'BurstOn';
EventName2 = 'PulseOn';
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName1);
epoch_pos2 = findcellstr(ST.events(:,1),EventName2);
if epoch_pos1 == 0 || epoch_pos2 == 0
    error('Event name not found');
end
stimes1 = ST.event_stimes{epoch_pos1};
stimes2 = ST.event_stimes{epoch_pos2};
time = win(1):dt:win(end);
valid_trials1 = find(~isnan(TE.(EventName1)));
% minfreq = min([TE.BurstNPulse]);
% maxpow = max([TE.PulsePower]);
% inx = ~isnan(TE.(EventName2)) & TE.BurstNPulse==minfreq & TE.PulsePower==maxpow;
% valid_trials2 = find(inx);
valid_trials2 = find(~isnan(TE.(EventName2)));
lm = 5000;    % downsaple if more pulses than 5000
if length(valid_trials2) > lm
    rp = randperm(length(valid_trials2));
    valid_trials2 = valid_trials2(sort(rp(1:lm)));
end

% Calculate bin rasters
spt1 = stimes2binraster(stimes1(valid_trials1),time,dt);
spt2 = stimes2binraster(stimes2(valid_trials2),time,dt);

% Set input arguments for rater plot and PSTH
if dsply
    SEvent = 'BurstOff';
    FNum = 2;
    parts = 'all';
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent = {{'BurstOff'}};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    
    % Plot raster plot and PSTH for 'BurstOn'
    figure
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName1,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
    
    % Plot raster plot and PSTH for 'PulseOn'
    figure
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName2,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% Calculate information distances and p values
res = 10;   % resolution in ms
dtt = dt * 1000;   % resolution of bin raster in ms
wn = win * 1000;   % window boundaries in ms
[p_value Idiff] = isikldist(spt1,spt2,dtt,wn,res);

% -------------------------------------------------------------------------
function [p_value Idiff] = isikldist(spt_baseline,spt_test,dt,win,res)

% Trial number and epoch length
[tno tl] = size(spt_baseline);

% Number of bins for ISI histograms
nmbn = round(res/dt);

% Pre-stimulus time window to consider for null hypothesis
st = abs(win(1)) / dt;   % number of pre-stim values in 'spt'

% ISI histogram - baseline
edges = 0:nmbn+1;
nm = floor(st/nmbn);
lsi = zeros(tno,nm);   % ISI's
slsi = zeros(tno,nm);  % sorted ISI's
hlsi = zeros(nmbn+1,nm);    % ISI hist.; +1: zero when no spike in the segment
nhlsi = zeros(nmbn+1,nm);   % normalized ISI histogram 
next = 1;
for t = 1:nmbn:st
    for k = 1:tno
        cspt = spt_baseline(k,t:t+nmbn-1);
        pki = find(cspt,1,'first');
        if ~isempty(pki)
            lsi(k,next) = pki;
        else
            lsi(k,next) = 0;
        end
    end
    slsi(:,next) = sort(lsi(:,next));
    hst = hist(slsi(:,next),edges);
    hlsi(:,next) = hst(1:end-1);
    nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));
    next = next + 1;
end

% ISI histogram - test
tno_test = size(spt_test,1);
lsi_tt = nan(tno_test,1);
for k = 1:tno_test
    cspt = spt_test(k,st+1:st+nmbn);
    pki = find(cspt,1,'first');
    if ~isempty(pki)
        lsi_tt(k,1) = pki;
    else
        lsi_tt(k,1) = 0;
    end
end
slsi_tt = sort(lsi_tt(:,1));
hst = hist(slsi_tt,edges);
hlsi(:,next) = hst(1:end-1);
nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));

% figure      % plot ISIs
% imagesc(lsi)
% figure      % plot sorted ISIs
% imagesc(slsi)
% figure      % plot ISI histograms
% imagesc(hlsi(2:end,:))

% Symmetric KL-divergence and JS-divergence
kn = st / nmbn + 1;
jsd = nan(kn,kn);  % pairwise modified JS-divergence (which is a metric!)
for k1 = 1:kn
    D1 = nhlsi(:,k1);
    for k2 = k1+1:kn
        D2 = nhlsi(:,k2);
        jsd(k1,k2) = sqrt(JSdiv(D1,D2)*2);
    end
end
% figure    % plot KL-distance
% imagesc(kld)

% Calculate p-value and information difference
[p_value Idiff] = makep(jsd,kn);
% keyboard

% -------------------------------------------------------------------------
function [p_value Idiff] = makep(kld,kn)
% Calculates p value from distance matrix.

pnhk = kld(1:kn-1,1:kn-1);
nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
testkld = median(kld(1:kn-1,kn));  % value to test
sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
p_value = length(find(nullhypkld>=testkld)) / sno;
Idiff = testkld - median(nullhypkld);