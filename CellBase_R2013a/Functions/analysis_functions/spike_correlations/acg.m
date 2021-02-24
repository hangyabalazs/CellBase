function [CCR, lags] = acg(cellids,window,varargin)
%ACG   Auto-correlation.
%   ACG(CELLIDS,WINDOW) calculates auto-correlations at 0.5 ms resolution 
%   for a given window size (WIN). For details on the algorithm, see XCORR. 
%
%   Required input arguments:
%       CELLIDS: cells included in the analysis.
%       WINDOW: ACG window size in seconds.
%
%   Optional input parameter-value pairs:
%       'dt', 0.0005 - ACG resolution (bin size)
%       'issave', false - controls saving behavior; plots and
%           auto-correlation matrices are only saved if 'issave' is 'true'
%       'resdir' - results directory. An 'ACG' directory is created and 
%           used on the analysis path if 'issave' is set to 'true' but no
%           results directory is specified
%           'nontetrodepairs' selects cells from other tetrodes,
%           'tetrodepairs' selects cells from the same tetrode and
%           'allpairs' selects all cells from the session
%       'minspikeno', 100 - minimal spike number to perform ACG calculation
%       'maxspikeno', 5000 - maximal spike number; the first N spikes are
%           included in the ACG calculation
%       'segfilter', 'none' - filter recording segments (see FINDSEGS3)
%           If set to 'none' (default), the entire recording is used
%       'filterinput', [] - additional input for FINDSEGS3
%       'longsegments', false - if true, only the longest segment is
%           analyzed after filtering
%       'seglim', 0.3 - if 'longsegments' is 'true', ACG is calculated if
%           the longest segment reaches 'seglim' (in seconds)
%
%   Examples for possible segment filters:
%         tseg = findSegs3(cell1,'segfilter','stimfb_excl_nb',...
%             'light_activation_duration',[-5 5],'margins',[0 0]);
%         tseg = findSegs3(cell1,'segfilter','prestim3');
%         tseg = findSegs3(cell1,'segfilter','fb_incl_nb',...
%             'feedback_duration',[-0.5 0.5],'margins',[0 0],'min_int',0);
%         tseg = findSegs3(cell1,'segfilter','cue_incl_nb',...
%             'feedback_duration',[-1.5 0],'margins',[0 0],'min_int',0);
%
%   Three indeces are derived from the auto-correlogram (ACG). Burst index
%   is calculated as the normalized difference between maximum ACG for lags
%   0-10 ms and mean ACG for lags 180-200 ms. The normalizing factor is the
%   greater of the two numbers, yielding and index between -1 and 1 (see
%   Royer et al., 2012). Theta index is calculated as the normalized
%   difference between mean ACG for a +-25 ms window around the peak within
%   lags 100 and 200 ms and the mean ACG for lags 180-200 and 65-85 ms.
%   Normalization is performed similarly as for the burst index. Refractory
%   is calculated as full width half hight of the central gap in the
%   smoothed ACG. Smoothing is performed by a 10 ms moving average. 
%
%   Reference:
%   Royer S, Zemelmen BV, Losonczy A, Kim J, Chance F, Magee JC, Buzsaki G
%   (2012) Control of timing, rate and bursts of hippocampal place cells by
%   dendritic and somatic inhibition. Nat Neurosci 15:769-775
%
%   See also CCG and XCORR.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addRequired(prs,'window',@isscalar)  % ACG window, in seconds
addParameter(prs,'dt',0.0005,@isnumeric)   % ACG resolution, in seconds
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'resdir',[],@(s)isdir(s)|isempty(s))   % results directory
addParameter(prs,'longsegments',false,@islogical)   % use only the longest segment after segment filtering
addParameter(prs,'seglim',0.3,@isnumeric);   % minimal segment length (s) to perform ACG calculation if 'longsegments' is 'true'
addParameter(prs,'segfilter','none',@(s)ischar(s)|iscellstr(s))   % filter segments
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'minspikeno',100,@isnumeric)   % calculate ACG above minimal spike number
addParameter(prs,'maxspikeno',5000,@isnumeric)   % restrict included spikes
parse(prs,cellids,window,varargin{:})
g = prs.Results;
if ischar(cellids)
    cellids = {cellids};  % one cell ID
end

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
if isempty(g.resdir)
    resdir = fullfile(DATAPATH,'ACG');  % results directory
    if ~isdir(resdir)
        mkdir(resdir)
    end
end
fnmm = 'ACG_matrices_.mat';   % filename for saving the result matrices

% Determine time window
sr = 1000;      % sampling rate
wn = g.window * sr;    % CCG window in data points
res = g.dt * sr;   % resolution for ACG in ms

% Cell loop for ACG
wb = waitbar(0,'Please wait...','Name','Running ACG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [g.minspikeno g.maxspikeno];   % include max 50000 spikes; calculate only if min 100 spikes by default
numCells = length(cellids);
[CCR, SCCR] = deal(zeros(numCells,2*wn/res));
[SegmentLength, BurstIndex, Refractory, ThetaIndex] = deal(nan(numCells,1));
for iC = 1:numCells   % loop through the cells
    cell = cellids{iC};
    if isequal(g.segfilter,'none')
        ncc = loadcb(cell,'SPIKES');   % use all spikes
    else
        tseg = findSegs3(cell,'segfilter',g.segfilter,...
            g.filterinput{:});  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if g.longsegments   % use the longest segment if it's longer than the threshold ('seglim')
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < g.seglim * sr
                continue
            end
        end
        SegmentLength(iC) = sum(ltseg);  % cumulative length of the segments
        ncc = extractSegSpikes(cell,tseg);   % find spikes in the time segments
    end
    
    % Implement upper spike number limits
    if length(ncc) > limit_spikes(2);      % crop if too long to avoid out of memory
        ncc = ncc(1:limit_spikes(2));
    end
    
    if length(ncc) > limit_spikes(1)     % implement minimum number of spikes
        [H1, ccr, lags] = acorr(ncc,wn,res);
        sccr = smooth(ccr,'linear',21);    % smoothed ACR
%         nqf = 1 / res * 1000 / 2;   % Nyquist freq.
%         flt = fir1(32,[4 10]/nqf,'bandpass');
%         sccr2 = filtfilt(flt,1,sccr);   % high-pass filter > 1 Hz
        hold on
        plot(lags,sccr,'Color',[0.7 0.7 0.7])
%         plot(lags,sccr2,'c')
        bar(lags(lags>-10&lags<10),ccr(lags>-10&lags<10),'FaceColor','g','EdgeColor','g')
        bar(lags(lags>-200&lags<-180),ccr(lags>-200&lags<-180),'FaceColor','r','EdgeColor','r')
        bar(lags(lags>180&lags<200),ccr(lags>180&lags<200),'FaceColor','r','EdgeColor','r')
        bar(lags((lags>=180&lags<=200)|(lags>=65&lags<=85)),...
            ccr((lags>=180&lags<=200)|(lags>=65&lags<=85)),...
            'FaceColor','c','EdgeColor','c');
        BurstIndex(iC) = burstinx(ccr,lags);   % burst index
        Refractory(iC) = refract(ccr,sccr,lags);   % refractory
        ThetaIndex(iC) = thetainx4(sccr,lags);   % theta index
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.95,regexprep(cell,'_',' '))
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,['Burst index: ' num2str(BurstIndex(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,['Refractory: ' num2str(Refractory(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.8,['Theta index: ' num2str(ThetaIndex(iC))])
        if g.issave   % save figure
            ncl = regexprep(cell,'\.','_');
            fnm = ['ACG_' ncl '.fig'];
            saveas(H1,fullfile(g.resdir,fnm))   % save CCG plot
            fnm2 = ['ACG_' ncl '.jpg'];
            saveas(H1,fullfile(g.resdir,fnm2))   % save CCG plot
            close(H1)
        end
        CCR(iC,:) = ccr;   % auto-correlogram
        SCCR(iC,:) = sccr;   % smoothed auto-correlogram
        disp(['Cell #' num2str(iC) ' / ' num2str(numCells) ' done......'])
    end
    
    % Save
    if g.issave
        if isequal(mod(iC,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(g.resdir,fnm),'cellids','CCR','SCCR','lags',...
                'SegmentLength','BurstIndex','Refractory','ThetaIndex')
            disp('Autosave done.')
        end
    end
    
    waitbar(iC/numCells)
end
close(wb)   % eliminate progress indicator

% Save
if g.issave
    save(fullfile(g.resdir,fnmm),'cellids','CCR','SCCR','lags',...
        'SegmentLength','BurstIndex','Refractory','ThetaIndex')
end

% -------------------------------------------------------------------------
function [H1, ccr, lags] = acorr(ncc,wn,res)

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;
mn = nc(1);  % only relative spike times count; avoid out of memory
nc = nc - mn;
nc(nc<res/2) = [];  % drop spikes in the first few ms - their indices would be rounded to zero below
wn2 = wn / 1000;    % window size in seconds

% Auto-correlogram
zunit = zeros(1,round(nc(end)/res)+5);
zunit(round(nc/res)) = 1;
[ccr, lags] = xcorr(zunit,zunit,wn2*sr/res);     % 1->2; window: -wn ms - wn ms
ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
lags(length(lags)/2+0.5) = [];
lags = lags * res;   % in ms

% Plot
H1 = figure;
bar(lags,ccr,'FaceColor','black')
set(gca,'XLim',[-wn wn])

% -------------------------------------------------------------------------
function bi = burstinx(ccr,lags)

a = max(ccr(lags>0&lags<10));
b = mean(ccr(lags>180&lags<200));
bi = (a - b) / max(a,b);   % burst index

% -------------------------------------------------------------------------
function tau = refract(ccr,sccr,lags)

mx = max(sccr);   % peak: from smoothed ACG
halfmx = mx / 2;   % half-max
ccr2 = ccr(lags>0);   % half ACG
lags2 = lags(lags>0);   % corresponding lag values
tau = lags2(find(ccr2>halfmx,1,'first'));   % refractory

% -------------------------------------------------------------------------
function thinx = thetainx4(ccr,lags)

% Theta index
thpeak = ccr(lags>=100&lags<=200);   % ACG first theta peak
[~, pl] = max(thpeak);
ploc = lags(find(lags>=100,1,'first')+pl-1);   % peak location
mthp = mean(ccr(lags>=ploc-25&lags<ploc+25));

thtrough = ccr((lags>=225&lags<=275)|(lags>=65&lags<=85));   % ACG first theta trough
mtht = mean(thtrough);
thinx = (mthp - mtht) / max(mthp,mtht);   % theta index