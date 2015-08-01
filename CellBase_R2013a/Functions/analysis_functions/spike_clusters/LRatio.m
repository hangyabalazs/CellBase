function [DM LrC valid_channels] = LRatio(cellid,feature_names,valid_channels)
%LRATIO   L-ratio and Isolation Distance.
%   [DM LRC] = LRATIO(CELLID,FEATURE_NAMES,VALID_CHANNELS) calculates
%   cluster quality measures (Isolation Distance, DM; L-ratio,LRC) for a
%   given cell (CELLID) using the given clustering features
%   (FEATURE_NAMES). The third input argument (VALID_CHANNELS) is optional;
%   it is a 4-element 0-1 array determining which of the four tetrode
%   channels were functional. If it is omitted, the program defines it
%   based on valid channels having nonzero waveform energy.
%
%   [DM LRC VALID_CHANNELS] = LRATIO(CELLID,FEATURE_NAMES,VALID_CHANNELS)
%   also returns channel validity.
%
%   See also MAHAL.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com

% Channel validity
if nargin < 3
    valid_channels = check_channel_validity(cellid);
end

% Parse cellID
[r,s,t,u] = cellid2tags(cellid);

% Load Ntt file
Nttfn = cellid2fnames(cellid,'Ntt');
all_spikes = LoadTT_NeuralynxNT(Nttfn);
TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
all_spikes = all_spikes * TIMEFACTOR;
spk = loadcb(cellid,'Spikes');
n = length(all_spikes);   % avoid subsampling of spikes due to rounding errors
spk = round(spk*1e6)/1e6;
all_spikes = round(all_spikes*1e6)/1e6;
[jnk, inx] = intersect(all_spikes,spk);
if ~isequal(jnk,spk)   % internal check for spike times
    error('LRatio:SpikeTimeMismatch','Mismatch between saved spike times and Ntt time stamps.')
end
cinx = setdiff(1:n,inx);

% Feature matrix
X = [];
for k = 1:length(feature_names)
    basename = [getpref('cellbase','cell_pattern') num2str(t)];
    propfn = [basename '_' feature_names{k}];   % name of feature file (e.g. TT1_Amplitude)
    sessionpath = cellid2fnames(cellid,'sess');
    propfn_path = [sessionpath filesep 'FD'];   % where the feature file can be found
    if ~isdir(propfn_path)
        propfn_path = sessionpath;
    end
    propfn_full = [propfn_path filesep propfn];   % full path of feature file
    try
        wf_prop = load([propfn_full '.fd'],'-mat');     % load feature file
    catch    %#ok<CTCH>
        disp('Calculating missing feature data.')       % calculate feature file if it was not found
        calculate_features(sessionpath,propfn_path,feature_names(k),basename,valid_channels)
        wf_prop = load([propfn_full '.fd'],'-mat');     % load feature file
    end
    wf_prop = wf_prop.FeatureData;
    
    if ~isequal(size(wf_prop,2),sum(valid_channels))
        wf_prop  = wf_prop(:,valid_channels);   % estimated and original valid_channels don't match
    end
    X = [X; wf_prop']; %#ok<AGROW>
end

% Isolation Distance
XC = X(:,inx);
D = mahal(X',XC');
DC = D(cinx);
sDC = sort(DC);
linx = length(inx);
if linx <= length(sDC)
    DM = sDC(length(inx));
else
    DM = Inf;   % more spikes in the cluster than outside of it
end

% L-ratio
if ~isempty(DC)
    df = size(X,1);
    LC = sum(1-chi2cdf(DC,df));
    nc = size(XC,2);
    LrC = LC / nc;
else
    LrC = Inf;  % no spikes outside of the cluster: multiunit
end

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